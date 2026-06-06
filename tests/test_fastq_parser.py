import logging

from modules.fastq_parser import (
    FastqRun,
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

    assert grouped["SAMPLE_A"] == [FastqRun("SAMPLE_A", "RUN_1", [r1, r2], None)]
    assert grouped["SAMPLE_B"] == [FastqRun("SAMPLE_B", "ERR123", [single], None)]


def test_group_fastqs_by_run_uses_lane_as_run_for_flat_illumina_files(tmp_path):
    r1_l1 = touch(tmp_path / "S1_L001_R1_001.fastq.gz")
    r2_l1 = touch(tmp_path / "S1_L001_R2_001.fastq.gz")
    r1_l2 = touch(tmp_path / "S1_L002_R1_001.fastq.gz")
    r2_l2 = touch(tmp_path / "S1_L002_R2_001.fastq.gz")

    grouped = group_fastqs_by_run(tmp_path)

    assert grouped["S1"] == [
        FastqRun("S1", "S1_L001", [r1_l1, r2_l1], None),
        FastqRun("S1", "S1_L002", [r1_l2, r2_l2], None),
    ]


def test_group_fastqs_by_run_supports_run_folder_with_multiple_samples(tmp_path):
    run1 = "210416_A01210_0064_BHWKFYDMXX"
    run2 = "210419_A00581_0145_AHWKL3DMXX"
    s1_run1_r1 = touch(tmp_path / run1 / "PE-AncientHorses-01_S1_L001_R1_001.fastq.gz")
    s1_run1_r2 = touch(tmp_path / run1 / "PE-AncientHorses-01_S1_L001_R2_001.fastq.gz")
    s1_run2_r1 = touch(tmp_path / run2 / "PE-AncientHorses-01_S1_L002_R1_001.fastq.gz")
    s1_run2_r2 = touch(tmp_path / run2 / "PE-AncientHorses-01_S1_L002_R2_001.fastq.gz")
    s2_run1_r1 = touch(tmp_path / run1 / "PE-AncientHorses-02_S2_L001_R1_001.fastq.gz")
    s2_run1_r2 = touch(tmp_path / run1 / "PE-AncientHorses-02_S2_L001_R2_001.fastq.gz")

    grouped = group_fastqs_by_run(tmp_path)

    assert grouped["PE-AncientHorses-01_S1"] == [
        FastqRun("PE-AncientHorses-01_S1", f"{run1}_L001", [s1_run1_r1, s1_run1_r2], None),
        FastqRun("PE-AncientHorses-01_S1", f"{run2}_L002", [s1_run2_r1, s1_run2_r2], None),
    ]
    assert grouped["PE-AncientHorses-02_S2"] == [
        FastqRun("PE-AncientHorses-02_S2", f"{run1}_L001", [s2_run1_r1, s2_run1_r2], None),
    ]


def test_group_fastqs_by_run_keeps_sample_folder_when_parent_matches_sample_hint(tmp_path):
    r1 = touch(tmp_path / "S1" / "S1_L001_R1_001.fastq.gz")
    r2 = touch(tmp_path / "S1" / "S1_L001_R2_001.fastq.gz")

    grouped = group_fastqs_by_run(tmp_path)

    assert grouped["S1"] == [FastqRun("S1", "S1_L001", [r1, r2], None)]


def test_group_fastqs_by_run_prefers_paired_when_single_is_mixed(tmp_path, caplog):
    r1 = touch(tmp_path / "S1_R1.fastq.gz")
    r2 = touch(tmp_path / "S1_R2.fastq.gz")
    touch(tmp_path / "S1.fastq.gz")

    with caplog.at_level(logging.WARNING):
        grouped = group_fastqs_by_run(tmp_path)

    assert grouped["S1"] == [FastqRun("S1", "S1", [r1, r2], None)]
    assert "single を無視" in caplog.text


def test_group_fastqs_by_run_keeps_multiple_sample_single_fastqs(tmp_path):
    fastq_a = touch(tmp_path / "SAMEA103910521" / "BER12_M__BER12_M_SCY1.1_user_GACGAC_EAV4_20141017_hiseq3a.fastq.gz")
    fastq_b = touch(tmp_path / "SAMEA103910521" / "BER12_M__BER12_M_SCY1.2_user_TCTCGC_EAV4_20141030_hiseq3b.fastq.gz")

    grouped = group_fastqs_by_run(tmp_path)

    assert grouped["SAMEA103910521"] == [
        FastqRun(
            "SAMEA103910521",
            "BER12_M__BER12_M_SCY1.1_user_GACGAC_EAV4_20141017_hiseq3a",
            [fastq_a],
            "SCY1.1",
        ),
        FastqRun(
            "SAMEA103910521",
            "BER12_M__BER12_M_SCY1.2_user_TCTCGC_EAV4_20141030_hiseq3b",
            [fastq_b],
            "SCY1.2",
        ),
    ]


def test_group_fastqs_by_run_pairs_multiple_reads_without_dropping(tmp_path):
    r1a = touch(tmp_path / "S1" / "RUN1" / "laneA_R1.fastq.gz")
    r2a = touch(tmp_path / "S1" / "RUN1" / "laneA_R2.fastq.gz")
    r1b = touch(tmp_path / "S1" / "RUN1" / "laneB_R1.fastq.gz")
    r2b = touch(tmp_path / "S1" / "RUN1" / "laneB_R2.fastq.gz")

    grouped = group_fastqs_by_run(tmp_path)

    assert grouped["S1"] == [
        FastqRun("S1", "RUN1", [r1a, r2a], None),
        FastqRun("S1", "RUN1_2", [r1b, r2b], None),
    ]


def test_group_fastqs_by_run_excludes_unpaired_reads(tmp_path, caplog):
    r1 = touch(tmp_path / "S1" / "RUN1" / "reads_R1.fastq.gz")

    with caplog.at_level(logging.WARNING):
        grouped = group_fastqs_by_run(tmp_path)

    assert "S1" not in grouped
    assert "片側だけの paired-end FASTQ を除外" in caplog.text


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
