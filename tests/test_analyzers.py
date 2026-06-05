import subprocess
from types import SimpleNamespace

import pytest

import modules.analyzers as analyzers
from modules.analyzers import HaplotypeCaller, MapDamageAnalyzer, QualimapAnalyzer
from tests.conftest import touch


def test_mapdamage_filters_indexes_and_runs_mapdamage(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    bam = tmp_path / "input.bam"
    commands = []

    def fake_run(cmd, **kwargs):
        commands.append((cmd, kwargs))
        return SimpleNamespace(stdout="")

    stream_calls = []
    monkeypatch.setattr(analyzers.subprocess, "run", fake_run)
    monkeypatch.setattr(analyzers, "_stream_cmd", lambda cmd, step, env=None: stream_calls.append((cmd, step, env)))

    out = MapDamageAnalyzer(cfg).run_mapdamage("S1", bam)

    assert out == cfg.results_dir / "S1" / "mapdamage"
    assert commands[0][1]["shell"] is True
    assert "samtools view -h -q 30" in commands[0][0]
    assert commands[1][0] == ["samtools", "index", str(cfg.temp_dir / "S1_filtered.sorted.bam")]
    assert stream_calls[0][0][:2] == ["mapDamage", "-i"]
    assert "--merge-libraries" in stream_calls[0][0]


def test_mapdamage_returns_none_when_filter_fails(monkeypatch, make_config, tmp_path):
    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, args[0])

    monkeypatch.setattr(analyzers.subprocess, "run", fail)

    assert MapDamageAnalyzer(make_config()).run_mapdamage("S1", tmp_path / "input.bam") is None


def test_qualimap_skips_missing_empty_and_zero_mapped_bam(monkeypatch, make_config, tmp_path):
    analyzer = QualimapAnalyzer(make_config())
    assert analyzer.run_qualimap("S1", tmp_path / "missing.bam") is None

    empty = tmp_path / "empty.bam"
    empty.write_bytes(b"")
    assert analyzer.run_qualimap("S1", empty) is None

    bam = touch(tmp_path / "reads.bam")
    monkeypatch.setattr(
        analyzers.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout="chr1\t100\t0\t0\n"),
    )
    assert analyzer.run_qualimap("S1", bam) is None


def test_qualimap_runs_with_java_tool_options(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    bam = touch(tmp_path / "reads.bam")
    monkeypatch.setattr(
        analyzers.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(stdout="chr1\t100\t3\t0\n"),
    )
    stream_calls = []
    monkeypatch.setattr(analyzers, "_stream_cmd", lambda cmd, step, env=None: stream_calls.append((cmd, step, env)))

    out = QualimapAnalyzer(cfg).run_qualimap("S1", bam)

    assert out == cfg.results_dir / "S1" / "qualimap"
    assert stream_calls[0][0][:2] == ["qualimap", "bamqc"]
    assert "-XX:+IgnoreUnrecognizedVMOptions" in stream_calls[0][2]["JAVA_TOOL_OPTIONS"]


def test_haplotypecaller_builds_gatk_command(monkeypatch, make_config, tmp_path):
    cfg = make_config()
    bam = tmp_path / "dedup.bam"
    stream_calls = []
    monkeypatch.setattr(analyzers, "_stream_cmd", lambda cmd, step, env=None: stream_calls.append((cmd, step, env)))

    out = HaplotypeCaller(cfg).run_haplotypecaller("S1", bam)

    assert out == cfg.results_dir / "S1" / "vcf_files" / "S1.vcf"
    cmd = stream_calls[0][0]
    assert cmd[:2] == ["gatk", "HaplotypeCaller"]
    assert cmd[cmd.index("-R") + 1] == str(cfg.reference_genome)
    assert cmd[cmd.index("-I") + 1] == str(bam)
    assert cmd[cmd.index("--native-pair-hmm-threads") + 1] == "4"


def test_haplotypecaller_returns_none_on_failure(monkeypatch, make_config, tmp_path):
    def fail(cmd, step, env=None):
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(analyzers, "_stream_cmd", fail)

    assert HaplotypeCaller(make_config()).run_haplotypecaller("S1", tmp_path / "dedup.bam") is None
