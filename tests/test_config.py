import argparse
import sys

import pytest

import config as config_module
from config import PipelineConfig, parse_args
from tests.conftest import touch


def make_args(tmp_path, **overrides):
    args = argparse.Namespace(
        base_dir=tmp_path / "base",
        fastq_dir=None,
        bam_dir=None,
        project_accession="PRJ",
        reference_genome=None,
        data_type="ancient",
        picard_jar=tmp_path / "picard.jar",
    )
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def test_pipeline_config_derives_project_from_input_dirs(tmp_path):
    bam_dir = tmp_path / "bams"
    bam_dir.mkdir()
    cfg = PipelineConfig(make_args(tmp_path, bam_dir=bam_dir))
    assert cfg.project_accession == "bams"

    fastq_dir = tmp_path / "fastqs"
    fastq_dir.mkdir()
    cfg = PipelineConfig(make_args(tmp_path, fastq_dir=fastq_dir))
    assert cfg.project_accession == "fastqs"


def test_pipeline_config_creates_expected_directories(tmp_path):
    cfg = PipelineConfig(make_args(tmp_path))

    assert cfg.raw_data_dir.is_dir()
    assert cfg.results_dir.is_dir()
    assert cfg.logs_dir.is_dir()
    assert cfg.temp_dir.is_dir()


def test_parse_args_rejects_exclusive_input_modes(monkeypatch, tmp_path):
    monkeypatch.setattr(
        sys,
        "argv",
        ["main.py", "--bam_dir", str(tmp_path / "bams"), "--fastq_dir", str(tmp_path / "fastqs")],
    )
    with pytest.raises(SystemExit):
        parse_args()

    monkeypatch.setattr(
        sys,
        "argv",
        ["main.py", "--bam_dir", str(tmp_path / "bams"), "--download-via-https"],
    )
    with pytest.raises(SystemExit):
        parse_args()


def test_parse_args_accepts_no_progress(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["main.py", "--no-progress"])

    assert parse_args().no_progress is True


def test_validate_environment_reports_missing_items(monkeypatch, tmp_path):
    ref = touch(tmp_path / "ref.fa")
    args = make_args(tmp_path, reference_genome=ref, picard_jar=tmp_path / "missing.jar")
    cfg = PipelineConfig(args)
    monkeypatch.setattr(config_module.shutil, "which", lambda tool: None)

    with pytest.raises(SystemExit):
        cfg.validate_environment()


def test_validate_environment_succeeds_when_tools_and_indexes_exist(monkeypatch, tmp_path):
    ref = touch(tmp_path / "ref.fa")
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        touch(tmp_path / ("ref.fa" + suffix))
    picard = touch(tmp_path / "picard.jar")
    cfg = PipelineConfig(make_args(tmp_path, reference_genome=ref, picard_jar=picard))
    monkeypatch.setattr(config_module.shutil, "which", lambda tool: f"/bin/{tool}")

    cfg.validate_environment()
