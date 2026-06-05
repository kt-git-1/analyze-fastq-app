from pathlib import Path

from modules.bwa_mapper import BWAMapper
from tests.conftest import touch


def test_adapter_removal_commands_for_ancient_and_modern(make_config, tmp_path):
    mapper = BWAMapper(make_config(data_type="ancient"))
    cmd = mapper._build_adapter_removal_cmd_paired("r1.fq.gz", "r2.fq.gz", tmp_path / "out")
    assert "--collapse" in cmd
    assert cmd[-5:] == ["--minquality", "20", "--minlength", "30", "--collapse"]

    mapper = BWAMapper(make_config(data_type="modern"))
    cmd = mapper._build_adapter_removal_cmd_single("r1.fq.gz", tmp_path / "out")
    assert "--file2" not in cmd
    assert "--collapse" not in cmd
    assert cmd[-4:] == ["--minquality", "25", "--minlength", "25"]


def test_bwa_pe_command_includes_threads_reference_and_read_group(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    mapper = BWAMapper(cfg)
    captured = {}
    monkeypatch.setattr(mapper, "_pipe_bwa_sort", lambda cmd, bam_out, label: captured.update(cmd=cmd, bam_out=bam_out, label=label) or bam_out)

    out = mapper._run_bwa_pe_and_sort(Path("r1.fq.gz"), Path("r2.fq.gz"), tmp_path / "out.bam", "@RG\\tID:run\\tSM:S1", "run")

    assert out == tmp_path / "out.bam"
    assert captured["cmd"][:2] == ["bwa", "mem"]
    assert captured["cmd"][captured["cmd"].index("-t") + 1] == "4"
    assert str(cfg.reference_genome) in captured["cmd"]
    assert "@RG\\tID:run\\tSM:S1" in captured["cmd"]
    assert captured["cmd"][-2:] == ["r1.fq.gz", "r2.fq.gz"]


def test_run_mapping_pipeline_modern_pe_requires_adapter_outputs(monkeypatch, make_config, tmp_path):
    cfg = make_config(data_type="modern")
    mapper = BWAMapper(cfg)
    monkeypatch.setattr(mapper, "_run_streaming", lambda cmd, step_name: None)

    assert mapper.run_mapping_pipeline("S1", "RUN1", [Path("r1"), Path("r2")]) is None


def test_run_mapping_pipeline_modern_pe_moves_outputs_and_maps(monkeypatch, make_config, tmp_path):
    cfg = make_config(data_type="modern")
    mapper = BWAMapper(cfg)
    monkeypatch.setattr(mapper, "_run_streaming", lambda cmd, step_name: None)
    temp_prefix = cfg.temp_dir / "S1" / "RUN1"
    touch(temp_prefix.with_suffix(".pair1.truncated"))
    touch(temp_prefix.with_suffix(".pair2.truncated"))
    monkeypatch.setattr(mapper, "_run_bwa_pe_and_sort", lambda r1, r2, out, rg, label: out)

    out = mapper.run_mapping_pipeline("S1", "RUN1", [Path("r1"), Path("r2")])

    assert out == cfg.results_dir / "S1" / "runs" / "RUN1" / "bam_files" / "RUN1.sorted.bam"


def test_ancient_pe_maps_available_outputs_and_merges(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    mapper = BWAMapper(cfg)
    prefix = tmp_path / "temp" / "RUN1"
    bam_dir = tmp_path / "bam"
    bam_dir.mkdir()
    touch(prefix.with_suffix(".collapsed.truncated"))
    touch(prefix.with_suffix(".pair1.truncated"))
    touch(prefix.with_suffix(".pair2.truncated"))
    calls = []
    monkeypatch.setattr(mapper, "_run_bwa_se_and_sort", lambda fastq, out, rg, label: calls.append(("se", fastq)) or touch(out))
    monkeypatch.setattr(mapper, "_run_bwa_pe_and_sort", lambda r1, r2, out, rg, label: calls.append(("pe", r1, r2)) or touch(out))
    monkeypatch.setattr(mapper, "_merge_bams", lambda bams, merged: calls.append(("merge", tuple(bams), merged)) or True)

    assert mapper._map_ancient_pe(prefix, bam_dir, "RUN1", "rg", "S1") == bam_dir / "RUN1.merged.sorted.bam"
    assert [call[0] for call in calls] == ["se", "pe", "merge"]
