from modules.bam_processor import BAMProcessor
from tests.conftest import touch


def test_process_run_bam_builds_cleansam_command(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    processor = BAMProcessor(cfg)
    commands = []
    monkeypatch.setattr(processor, "_run_cmd", lambda cmd, step: commands.append((cmd, step)))

    out = processor.process_run_bam("S1", "RUN1", tmp_path / "soft.bam")

    assert out == cfg.results_dir / "S1" / "runs" / "RUN1" / "bam_files" / "RUN1.clean.bam"
    assert commands[0][1] == "Picard CleanSam"
    assert "CleanSam" in commands[0][0]
    assert "I=" + str(tmp_path / "soft.bam") in commands[0][0]


def test_merge_and_dedup_builds_commands_for_multiple_bams(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    processor = BAMProcessor(cfg)
    bam1 = touch(tmp_path / "run1.bam")
    bam2 = touch(tmp_path / "run2.bam")
    commands = []
    monkeypatch.setattr(processor, "_run_cmd", lambda cmd, step: commands.append((cmd, step)))

    out = processor.merge_and_dedup("S1", [bam1, bam2])

    assert out == cfg.results_dir / "S1" / "dedup" / "S1.dedup.sorted.bam"
    assert commands[0][0][:3] == ["samtools", "merge", "-f"]
    assert commands[1][1] == "Picard MarkDuplicates"
    assert commands[2][0][:3] == ["samtools", "sort", "-o"]
    assert commands[3][0][:2] == ["samtools", "index"]


def test_merge_and_dedup_returns_none_without_valid_bams(make_config):
    assert BAMProcessor(make_config()).merge_and_dedup("S1", []) is None
