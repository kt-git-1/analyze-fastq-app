import pytest

from modules.bam_parser import group_bams_by_sample, infer_sample_name_from_bam
from tests.conftest import touch


@pytest.mark.parametrize(
    ("filename", "expected"),
    [
        ("SAMPLE_A.dedup.sorted.bam", "SAMPLE_A"),
        ("SAMPLE_B.sorted.bam", "SAMPLE_B"),
        ("SAMPLE_C.dedup.bam", "SAMPLE_C"),
        ("SAMPLE_D.bam", "SAMPLE_D"),
    ],
)
def test_infer_sample_name_from_bam(filename, expected):
    assert infer_sample_name_from_bam(filename) == expected


def test_group_bams_by_sample_uses_pattern_and_recurses(tmp_path):
    bam = touch(tmp_path / "nested" / "S1.dedup.sorted.bam")
    touch(tmp_path / "nested" / "S2.raw.bam")

    assert group_bams_by_sample(tmp_path, "*.dedup.sorted.bam") == {"S1": bam}


def test_group_bams_by_sample_rejects_duplicate_inferred_sample(tmp_path):
    touch(tmp_path / "S1.bam")
    touch(tmp_path / "nested" / "S1.dedup.sorted.bam")

    with pytest.raises(ValueError, match="同一サンプル名"):
        group_bams_by_sample(tmp_path)
