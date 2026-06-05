import argparse
import hashlib
import sys
import threading
import time
import types
from pathlib import Path
from types import SimpleNamespace
from unittest import mock
import unittest

# The production pipeline depends on requests/tqdm at import time, but these
# tests mock network/progress behavior and should run in a minimal environment.
if "requests" not in sys.modules:
    requests_stub = types.ModuleType("requests")
    requests_stub.Session = lambda: object()
    requests_stub.exceptions = types.SimpleNamespace(RequestException=Exception)
    sys.modules["requests"] = requests_stub
if "pysam" not in sys.modules:
    sys.modules["pysam"] = types.ModuleType("pysam")
if "tqdm" not in sys.modules:
    tqdm_stub = types.ModuleType("tqdm")

    class _Tqdm:
        def __init__(self, iterable=None, *args, **kwargs):
            self.iterable = iterable

        def __iter__(self):
            return iter(self.iterable or [])

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def update(self, n=1):
            return None

    tqdm_stub.tqdm = _Tqdm
    sys.modules["tqdm"] = tqdm_stub

import main
from config import PipelineConfig, parse_args
from modules.fastq_parser import group_fastqs_by_run


def _args(tmp_path: Path, **overrides):
    values = dict(
        project_accession="PRJEB19970",
        base_dir=tmp_path,
        reference_genome=tmp_path / "reference.fa",
        workers=2,
        download_protocol="http",
        max_retries=1,
        threads=1,
        java_mem="1g",
        fastq_dir=None,
        data_type="ancient",
        force=False,
        picard_jar=tmp_path / "picard.jar",
        rg_library="unknown",
        rg_center="unknown",
        parallel_samples=1,
        download_via_https=False,
        download_only=False,
    )
    values.update(overrides)
    return argparse.Namespace(**values)


class PipelineFeatureTests(unittest.TestCase):
    def test_download_only_cli_option_is_available(self):
        with mock.patch("sys.argv", ["main.py", "--project_accession", "PRJEB19970", "--download-only"]):
            args = parse_args()

        self.assertEqual(args.project_accession, "PRJEB19970")
        self.assertTrue(args.download_only)

    def test_main_download_only_returns_before_analysis_validation(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as td:
            tmp_path = Path(td)
            args = _args(tmp_path, download_only=True)

            with mock.patch("main.parse_args", return_value=args), \
                 mock.patch.object(PipelineConfig, "validate_environment") as validate_environment, \
                 mock.patch("main.download_project_fastqs", return_value=tmp_path / "downloads") as download_project_fastqs, \
                 mock.patch("main.collect_fastqs") as collect_fastqs, \
                 mock.patch("main.run_sample_analyses") as run_sample_analyses:
                main.main()

            download_project_fastqs.assert_called_once()
            validate_environment.assert_not_called()
            collect_fastqs.assert_not_called()
            run_sample_analyses.assert_not_called()

    def test_main_analysis_mode_runs_fastq_collection_and_sample_analysis(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as td:
            tmp_path = Path(td)
            reference = tmp_path / "reference.fa"
            reference.write_text(">chr1\nACGT\n")
            reference.with_suffix(".fai").write_text("chr1\t4\t6\t4\t5\n")
            args = _args(tmp_path, reference_genome=reference, parallel_samples=2)
            sample_to_runs = {"SAMEA1": [("ERR1", [tmp_path / "SAMEA1_R1.fastq.gz"])]}

            with mock.patch("main.parse_args", return_value=args), \
                 mock.patch.object(PipelineConfig, "validate_environment") as validate_environment, \
                 mock.patch("main.collect_fastqs", return_value=sample_to_runs) as collect_fastqs, \
                 mock.patch("main.run_sample_analyses", return_value=(["SAMEA1"], [])) as run_sample_analyses:
                main.main()

            validate_environment.assert_called_once()
            collect_fastqs.assert_called_once()
            run_sample_analyses.assert_called_once()
            self.assertEqual(run_sample_analyses.call_args.args[2], 2)

    def test_download_project_fastqs_uses_prjeb19970_and_writes_ena_fastq_tree(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as td:
            tmp_path = Path(td)
            config = PipelineConfig(_args(tmp_path))
            payload = b"@r1\nACGT\n+\n!!!!\n"
            md5 = hashlib.md5(payload).hexdigest()
            created_raw = []

            class FakeDownloader:
                def __init__(self, cfg):
                    self.config = cfg

                def get_api_response(self, project_accession, session):
                    self.project_accession = project_accession
                    self.session = session
                    self.config.seen_project_accession = project_accession
                    return "sample_accession\tsubmitted_ftp\tsubmitted_md5\n"

                def parse_response_with_checksums(self, response_data):
                    self.config.seen_response = response_data
                    return {"SAMEA000001": [("ftp.sra.ebi.ac.uk/vol1/fastq/ERR000/ERR000001/test.fastq.gz", md5)]}

                def download_sample_data(self, sample_acc, ftp_urls):
                    raw = self.config.raw_data_dir / sample_acc / "test.fastq.gz"
                    raw.parent.mkdir(parents=True, exist_ok=True)
                    raw.write_bytes(payload)
                    created_raw.append(raw)
                    return [raw]

            with mock.patch("main.ENADownloader", FakeDownloader):
                out_dir = main.download_project_fastqs(config, session=object())

            self.assertEqual(config.seen_project_accession, "PRJEB19970")
            self.assertEqual(out_dir, config.results_dir / "ena_fastq")
            self.assertFalse(created_raw[0].exists())
            downloaded = out_dir / "SAMEA000001" / "test.fastq.gz"
            self.assertTrue(downloaded.exists())
            self.assertEqual(downloaded.read_bytes(), payload)

    def test_parallel_sample_analysis_runs_samples_concurrently(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as td:
            tmp_path = Path(td)
            config = SimpleNamespace(results_dir=tmp_path, args=SimpleNamespace(force=False))
            sample_to_runs = {
                "SAMEA1": [("ERR1", [tmp_path / "a.fastq.gz"])],
                "SAMEA2": [("ERR2", [tmp_path / "b.fastq.gz"])],
            }
            lock = threading.Lock()
            active = 0
            max_active = 0

            def fake_process_sample(sample_acc, runs, cfg):
                nonlocal active, max_active
                with lock:
                    active += 1
                    max_active = max(max_active, active)
                time.sleep(0.05)
                with lock:
                    active -= 1
                return sample_acc, True, ""

            with mock.patch("main.process_sample", side_effect=fake_process_sample):
                succeeded, failed = main.run_sample_analyses(sample_to_runs, config, parallel_samples=2)

            self.assertEqual(set(succeeded), {"SAMEA1", "SAMEA2"})
            self.assertEqual(failed, [])
            self.assertGreaterEqual(max_active, 2)

    def test_fastq_auto_detection_groups_paired_single_and_lane_runs(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as td:
            root = Path(td)
            for name in [
                "SAMPLE_A_L001_R1_001.fastq.gz",
                "SAMPLE_A_L001_R2_001.fastq.gz",
                "SAMPLE_A_L002_R1_001.fastq.gz",
                "SAMPLE_A_L002_R2_001.fastq.gz",
                "SAMPLE_B_R1.fastq.gz",
                "SAMPLE_B_R2.fastq.gz",
                "SAMPLE_C.fastq.gz",
            ]:
                (root / name).write_text("@r\nA\n+\n!\n")

            grouped = group_fastqs_by_run(root)

            self.assertEqual(len(grouped["SAMPLE_A"]), 2)
            sample_a_runs = {run_id: [p.name for p in paths] for run_id, paths in grouped["SAMPLE_A"]}
            self.assertEqual(sample_a_runs["SAMPLE_A_L001"], [
                "SAMPLE_A_L001_R1_001.fastq.gz",
                "SAMPLE_A_L001_R2_001.fastq.gz",
            ])
            self.assertEqual(sample_a_runs["SAMPLE_A_L002"], [
                "SAMPLE_A_L002_R1_001.fastq.gz",
                "SAMPLE_A_L002_R2_001.fastq.gz",
            ])
            self.assertEqual([p.name for _, paths in grouped["SAMPLE_B"] for p in paths], [
                "SAMPLE_B_R1.fastq.gz",
                "SAMPLE_B_R2.fastq.gz",
            ])
            self.assertEqual([p.name for _, paths in grouped["SAMPLE_C"] for p in paths], [
                "SAMPLE_C.fastq.gz",
            ])

    def test_completed_sample_is_skipped_without_instantiating_analysis_steps(self):
        from tempfile import TemporaryDirectory

        with TemporaryDirectory() as td:
            tmp_path = Path(td)
            config = SimpleNamespace(results_dir=tmp_path, args=SimpleNamespace(force=False))
            done_flag = tmp_path / "SAMEA_DONE" / ".done"
            done_flag.parent.mkdir(parents=True)
            done_flag.touch()

            with mock.patch("main.BWAMapper", side_effect=AssertionError("analysis should be skipped")):
                sample_acc, ok, step = main.process_sample(
                    "SAMEA_DONE",
                    [("ERR_DONE", [tmp_path / "done.fastq.gz"])],
                    config,
                )

            self.assertEqual(sample_acc, "SAMEA_DONE")
            self.assertTrue(ok)
            self.assertEqual(step, "")


if __name__ == "__main__":
    unittest.main()
