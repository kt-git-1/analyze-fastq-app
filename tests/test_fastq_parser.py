import logging

from modules.fastq_parser import (
    _classify_fastq_name,
    group_fastqs_by_run,
    merge_lanes_by_cat,
)
from tests.conftest import touch


def test_classify_supported_fastq_names():
    assert _classify_fastq_name("SAMPLE_A_L001_R1_001.fastq.gz") == (
        "SAMPLE_A",
        "R1",
        "SAMPLE_A_L001",
    )
    assert _classify_fastq_name("Horse01_L1_2.fq.gz") == ("Horse01", "R2", "Horse01_L1")
    assert _classify_fastq_name("sampleB-R1.fastq") == ("sampleB", "R1", None)
    assert _classify_fastq_name("sampleC_2.fq") == ("sampleC", "R2", None)
    assert _classify_fastq_name("ERR123456.fastq.gz") == (None, "single", None)


def test_group_fastqs_by_run_prefers_directory_sample_and_run(tmp_path):
    r1 = touch(tmp_path / "SAMPLE_A" / "RUN_1" / "reads_R1.fastq.gz")
    r2 = touch(tmp_path / "SAMPLE_A" / "RUN_1" / "reads_R2.fastq.gz")
    single = touch(tmp_path / "SAMPLE_B" / "ERR123.fastq.gz")

    grouped = group_fastqs_by_run(tmp_path)

    assert grouped["SAMPLE_A"] == [("RUN_1", [r1, r2])]
    assert grouped["SAMPLE_B"] == [("SAMPLE_B", [single])]


def test_group_fastqs_by_run_uses_lane_as_run_for_flat_illumina_files(tmp_path):
    r1_l1 = touch(tmp_path / "S1_L001_R1_001.fastq.gz")
    r2_l1 = touch(tmp_path / "S1_L001_R2_001.fastq.gz")
    r1_l2 = touch(tmp_path / "S1_L002_R1_001.fastq.gz")
    r2_l2 = touch(tmp_path / "S1_L002_R2_001.fastq.gz")

    grouped = group_fastqs_by_run(tmp_path)

    assert grouped["S1"] == [
        ("S1_L001", [r1_l1, r2_l1]),
        ("S1_L002", [r1_l2, r2_l2]),
    ]


def test_group_fastqs_by_run_prefers_paired_when_single_is_mixed(tmp_path, caplog):
    r1 = touch(tmp_path / "S1_R1.fastq.gz")
    r2 = touch(tmp_path / "S1_R2.fastq.gz")
    touch(tmp_path / "S1.fastq.gz")

    with caplog.at_level(logging.WARNING):
        grouped = group_fastqs_by_run(tmp_path)

    assert grouped["S1"] == [("S1", [r1, r2])]
    assert "single を無視" in caplog.text


def test_merge_lanes_by_cat_concatenates_pairs_and_treats_lone_pair_as_single(tmp_path):
    r1a = touch(tmp_path / "in" / "S1_L001_R1.fastq.gz", b"r1a")
    r1b = touch(tmp_path / "in" / "S1_L002_R1.fastq.gz", b"r1b")
    r2a = touch(tmp_path / "in" / "S1_L001_R2.fastq.gz", b"r2a")
    lone = touch(tmp_path / "in" / "S2_R1.fastq.gz", b"lone")
    reads = {
        "S1": {"R1": [r1a, r1b], "R2": [r2a], "single": []},
        "S2": {"R1": [lone], "R2": [], "single": []},
    }

    merged = merge_lanes_by_cat(reads, tmp_path / "merged", logging.getLogger(__name__))

    assert merged["S1"][0].read_bytes() == b"r1ar1b"
    assert merged["S1"][1] == r2a
    assert merged["S2"] == [lone]
