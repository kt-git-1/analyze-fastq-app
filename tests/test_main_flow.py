from pathlib import Path

import main as main_module
from modules.fastq_parser import FastqRun
from tests.conftest import touch


class FakeBWAMapper:
    result = Path("mapped.bam")
    calls = []

    def __init__(self, config):
        pass

    def run_mapping_pipeline(self, sample_acc, run_id, fastq_files, rg_library=None):
        self.calls.append((sample_acc, run_id, fastq_files, rg_library))
        return self.result


class FakeSoftClipper:
    result = Path("soft.bam")

    def __init__(self, config):
        pass

    def run_softclipping(self, sample_acc, run_id, bam_file):
        return self.result


class FakeBAMProcessor:
    run_result = Path("clean.bam")
    dedup_result = Path("dedup.bam")
    merge_calls = []

    def __init__(self, config):
        pass

    def process_run_bam(self, sample_acc, run_id, softclipped_bam):
        return self.run_result

    def merge_and_dedup(self, sample_acc, run_clean_bams):
        self.merge_calls.append((sample_acc, run_clean_bams))
        return self.dedup_result


class FakeMapDamageAnalyzer:
    result = Path("mapdamage")

    def __init__(self, config):
        pass

    def run_mapdamage(self, sample_acc, bam):
        return self.result


class FakeQualimapAnalyzer:
    result = Path("qualimap")

    def __init__(self, config):
        pass

    def run_qualimap(self, sample_acc, bam):
        return self.result


class FakeHaplotypeCaller:
    result = Path("sample.vcf")

    def __init__(self, config):
        pass

    def run_haplotypecaller(self, sample_acc, bam):
        return self.result


def patch_pipeline_classes(monkeypatch):
    monkeypatch.setattr(main_module, "BWAMapper", FakeBWAMapper)
    monkeypatch.setattr(main_module, "SoftClipper", FakeSoftClipper)
    monkeypatch.setattr(main_module, "BAMProcessor", FakeBAMProcessor)
    monkeypatch.setattr(main_module, "MapDamageAnalyzer", FakeMapDamageAnalyzer)
    monkeypatch.setattr(main_module, "QualimapAnalyzer", FakeQualimapAnalyzer)
    monkeypatch.setattr(main_module, "HaplotypeCaller", FakeHaplotypeCaller)
    monkeypatch.setattr(main_module, "cleanup_intermediate_file", lambda path, logger: None)


def reset_fake_results():
    FakeBWAMapper.result = Path("mapped.bam")
    FakeBWAMapper.calls = []
    FakeSoftClipper.result = Path("soft.bam")
    FakeBAMProcessor.run_result = Path("clean.bam")
    FakeBAMProcessor.dedup_result = Path("dedup.bam")
    FakeBAMProcessor.merge_calls = []
    FakeMapDamageAnalyzer.result = Path("mapdamage")
    FakeQualimapAnalyzer.result = Path("qualimap")
    FakeHaplotypeCaller.result = Path("sample.vcf")


def test_process_sample_skips_done_when_not_forced(make_config):
    cfg = make_config()
    (cfg.results_dir / "S1").mkdir(parents=True)
    (cfg.results_dir / "S1" / ".done").touch()

    assert main_module.process_sample("S1", [FastqRun("S1", "RUN1", [Path("r1")])], cfg) == ("S1", True, "")


def test_process_sample_success_creates_done(monkeypatch, make_config):
    reset_fake_results()
    patch_pipeline_classes(monkeypatch)
    cfg = make_config()

    assert main_module.process_sample("S1", [FastqRun("S1", "RUN1", [Path("r1")])], cfg) == ("S1", True, "")
    assert (cfg.results_dir / "S1" / ".done").exists()


def test_process_sample_maps_all_fastq_runs_and_merges_outputs(monkeypatch, make_config):
    reset_fake_results()
    patch_pipeline_classes(monkeypatch)
    cfg = make_config()
    runs = [
        FastqRun("S1", "RUN1", [Path("r1.fastq.gz")], "SCY1.1"),
        FastqRun("S1", "RUN2", [Path("r2.fastq.gz")], "SCY1.2"),
    ]

    assert main_module.process_sample("S1", runs, cfg) == ("S1", True, "")

    assert FakeBWAMapper.calls == [
        ("S1", "RUN1", [Path("r1.fastq.gz")], "SCY1.1"),
        ("S1", "RUN2", [Path("r2.fastq.gz")], "SCY1.2"),
    ]
    assert FakeBAMProcessor.merge_calls == [("S1", [Path("clean.bam"), Path("clean.bam")])]


def test_process_sample_reports_mapping_failure(monkeypatch, make_config):
    reset_fake_results()
    patch_pipeline_classes(monkeypatch)
    FakeBWAMapper.result = None

    assert main_module.process_sample("S1", [FastqRun("S1", "RUN1", [Path("r1")])], make_config()) == (
        "S1",
        False,
        "all runs failed",
    )


def test_process_sample_reports_downstream_failure(monkeypatch, make_config):
    reset_fake_results()
    patch_pipeline_classes(monkeypatch)
    FakeQualimapAnalyzer.result = None

    assert main_module.process_sample("S1", [FastqRun("S1", "RUN1", [Path("r1")])], make_config()) == (
        "S1",
        False,
        "Qualimap",
    )


def test_process_sample_from_dedup_bam_checks_input_and_runs_analyzers(monkeypatch, make_config, tmp_path):
    reset_fake_results()
    patch_pipeline_classes(monkeypatch)
    cfg = make_config()
    bam = touch(tmp_path / "S1.dedup.sorted.bam")
    monkeypatch.setattr(main_module, "_ensure_bam_index", lambda path: True)

    assert main_module.process_sample_from_dedup_bam("S1", bam, cfg) == ("S1", True, "")
    assert (cfg.results_dir / "S1" / ".done").exists()


def test_process_sample_from_dedup_bam_reports_missing_empty_and_index_failure(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    assert main_module.process_sample_from_dedup_bam("S1", tmp_path / "missing.bam", cfg) == (
        "S1",
        False,
        "BAM missing",
    )

    empty = tmp_path / "empty.bam"
    empty.write_bytes(b"")
    assert main_module.process_sample_from_dedup_bam("S1", empty, cfg) == ("S1", False, "BAM empty")

    bam = touch(tmp_path / "S1.bam")
    monkeypatch.setattr(main_module, "_ensure_bam_index", lambda path: False)
    assert main_module.process_sample_from_dedup_bam("S1", bam, cfg) == ("S1", False, "samtools index")
