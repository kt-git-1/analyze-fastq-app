from pathlib import Path
from typing import Dict


_BAM_SUFFIXES = (
    ".dedup.sorted.bam",
    ".sorted.bam",
    ".dedup.bam",
    ".bam",
)


def infer_sample_name_from_bam(bam_path: Path) -> str:
    """Infer sample name from a BAM filename."""
    name = Path(bam_path).name
    for suffix in _BAM_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def group_bams_by_sample(bam_dir: Path, pattern: str = "*.bam") -> Dict[str, Path]:
    """
    Recursively find BAM files and map them to sample names.

    Raises
    ------
    ValueError
        If multiple BAM files resolve to the same sample name.
    """
    bam_dir = Path(bam_dir)
    sample_to_bam: Dict[str, Path] = {}

    for bam_path in sorted(bam_dir.rglob(pattern)):
        if not bam_path.is_file():
            continue

        sample_acc = infer_sample_name_from_bam(bam_path)
        if not sample_acc:
            raise ValueError(f"サンプル名を推定できません: {bam_path}")

        existing = sample_to_bam.get(sample_acc)
        if existing is not None:
            raise ValueError(
                f"同一サンプル名に複数の BAM が見つかりました: "
                f"{sample_acc} ({existing}, {bam_path})"
            )

        sample_to_bam[sample_acc] = bam_path

    return sample_to_bam
