from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def make_config(tmp_path):
    def _make_config(**overrides):
        args = SimpleNamespace(
            threads=4,
            workers=2,
            java_mem="2g",
            picard_jar=tmp_path / "picard.jar",
            rg_library="lib",
            rg_center="center",
            force=False,
            download_protocol="http",
            max_retries=2,
            bam_pattern="*.bam",
            no_progress=True,
        )
        config = SimpleNamespace(
            args=args,
            base_dir=tmp_path,
            raw_data_dir=tmp_path / "raw_data" / "PROJECT",
            results_dir=tmp_path / "results" / "PROJECT",
            logs_dir=tmp_path / "logs",
            temp_dir=tmp_path / "temp",
            reference_genome=tmp_path / "ref.fa",
            data_type="ancient",
            fastq_dir=None,
            bam_dir=None,
            project_accession="PROJECT",
        )
        for path in (config.raw_data_dir, config.results_dir, config.logs_dir, config.temp_dir):
            path.mkdir(parents=True, exist_ok=True)
        for key, value in overrides.items():
            setattr(config, key, value)
        return config

    return _make_config


def touch(path: Path, content: bytes = b"x") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path
