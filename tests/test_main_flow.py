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

    def run_softclipping(self, sample_acc, run_id, bam_file, progress_enabled=True):
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
    calls = []

    def __init__(self, config):
        pass

    def run_haplotypecaller(self, sample_acc, bam):
        self.calls.append((sample_acc, bam))
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
    FakeHaplotypeCaller.calls = []


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


def test_process_sample_reuses_existing_vcf(monkeypatch, make_config):
    reset_fake_results()
    patch_pipeline_classes(monkeypatch)
    cfg = make_config()
    existing_vcf = cfg.results_dir / "S1" / "vcf_files" / "S1.vcf"
    touch(existing_vcf)

    assert main_module.process_sample("S1", [FastqRun("S1", "RUN1", [Path("r1")])], cfg) == ("S1", True, "")
    assert FakeHaplotypeCaller.calls == []


def test_cleanup_completed_fastq_intermediates_keeps_final_dedup_bam(make_config):
    cfg = make_config()
    sample_dir = cfg.results_dir / "S1"
    runs_file = touch(sample_dir / "runs" / "RUN1" / "bam_files" / "RUN1.sorted.bam")
    merged_bam = touch(sample_dir / "dedup" / "S1.merged.bam")
    marked_bam = touch(sample_dir / "dedup" / "S1.marked.bam")
    final_bam = touch(sample_dir / "dedup" / "S1.dedup.sorted.bam")
    final_bai = touch(sample_dir / "dedup" / "S1.dedup.sorted.bam.bai")
    temp_file = touch(cfg.temp_dir / "S1" / "RUN1.pair1.truncated")

    main_module._cleanup_completed_fastq_intermediates(cfg, ["S1"])

    assert not runs_file.exists()
    assert not merged_bam.exists()
    assert not marked_bam.exists()
    assert not temp_file.exists()
    assert final_bam.exists()
    assert final_bai.exists()


def test_write_sample_qc_summary(make_config):
    cfg = make_config()
    sample_dir = cfg.results_dir / "S1"
    touch(sample_dir / ".done")
    touch(sample_dir / "dedup" / "S1.dedup.sorted.bam")
    touch(sample_dir / "dedup" / "S1.dedup.sorted.bam.bai")
    touch(sample_dir / "mapdamage" / "Runtime_log.txt")
    touch(
        sample_dir / "qualimap" / "genome_results.txt",
        b"""number of reads = 1,000
number of mapped reads = 250 (25%)
duplication rate = 3.43%
mean insert size = 109.1
median insert size = 95
mean mapping quality = 53.4
GC percentage = 42.38%
general error rate = 0.0074
mean coverageData = 1.0108X
There is a 57.1% of reference with a coverageData >= 1X
There is a 22.37% of reference with a coverageData >= 2X
There is a 7% of reference with a coverageData >= 3X
""",
    )
    touch(
        sample_dir / "vcf_files" / "S1.vcf",
        b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tG\t30\t.\t.\n",
    )

    out = main_module._write_sample_qc_summary(cfg, ["S1"])

    text = out.read_text()
    assert "sample\tdata_type\tdone" in text
    assert "S1\tancient\ttrue" in text
    assert "\t1\t1000\t250\t1.0108\t57.1\t22.37\t7\t53.4" in text
    assert "pseudo-haploid cohort PCA" in text


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


def test_analysis_dashboard_renders_human_readable_status(monkeypatch):
    monkeypatch.setattr(main_module.time, "monotonic", lambda: 100.0)
    dashboard = main_module.AnalysisDashboard(
        "PRJEB19970",
        "FASTQ",
        total_samples=14,
        parallel_samples=4,
        threads_per_sample=8,
        enabled=False,
    )
    dashboard.started_at = 80.0
    dashboard.start_sample("SAMEA103910511", 10, "FASTQ  2 runs")
    dashboard.start_step("SAMEA103910511", "BWA mapping (RUN1) [1/2]", 1, 10)
    dashboard.finish_sample("SAMEA103910511", True, "")

    text = dashboard.render_text()

    assert "解析パイプライン: PRJEB19970" in text
    assert "全体進捗" in text
    assert "1/14 samples" in text
    assert "サンプル  SAMEA103910511" in text
    assert "現在      BWA mapping (RUN1) [1/2]" in text
    assert "ステップ" in text
    assert "1/10" in text
    assert "入力      FASTQ  2 runs" in text
    assert "並列数    4 samples" in text
    assert "スレッド  8 / sample" in text
    assert "最近のイベント" in text
    assert "解析完了" in text
    assert "root:" not in text
    assert "__main__:" not in text


def test_dashboard_screen_line_count_accounts_for_wrapped_and_wide_text():
    text = "解析パイプライン: ancient\n" + ("A" * 25)

    assert main_module._screen_line_count(text, 10) == 6
