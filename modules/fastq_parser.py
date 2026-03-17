"""
FASTQ ファイルの検出・分類・レーンマージを行うユーティリティ。

主なエントリポイント:
  - parse_fastq_general(fastq_dir) : ディレクトリ内の FASTQ をサンプル/R1/R2/single に分類
  - merge_lanes_by_cat(sample_to_reads, merged_dir, logger) : 複数レーンを1サンプル1セットに統合
"""

import logging
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional


def _strip_fastq_suffix(filename: str) -> str:
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def _classify_fastq_name(filename: str) -> tuple[Optional[str], str]:
    """
    FASTQ ファイル名から sample 名のヒントと read 種別を推定します。
    """
    stem = _strip_fastq_suffix(filename)
    patterns = [
        re.compile(r"^(?P<sample>.+?)_L\d{3}_R(?P<read>[12])(?:_\d+)?$"),
        re.compile(r"^(?P<sample>.+?)_L\d+_(?P<read>[12])$"),
        re.compile(r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]?\d+)?$", re.IGNORECASE),
        re.compile(r"^(?P<sample>.+?)[._-](?P<read>[12])$"),
    ]

    for pattern in patterns:
        match = pattern.match(stem)
        if match:
            return match.group("sample"), f"R{match.group('read')}"

    return None, "single"


def parse_fastq_general(fastq_dir: Path) -> Dict[str, Dict[str, List[Path]]]:
    """
    FASTQ を再帰探索して sample ごとに R1/R2/single に分類します。
    """
    fastq_dir = Path(fastq_dir)
    sample_to_reads: Dict[str, Dict[str, List[Path]]] = {}

    fastq_files = sorted(
        [
            *fastq_dir.rglob("*.fastq"),
            *fastq_dir.rglob("*.fastq.gz"),
            *fastq_dir.rglob("*.fq"),
            *fastq_dir.rglob("*.fq.gz"),
        ]
    )

    for fastq_file in fastq_files:
        rel_parts = fastq_file.relative_to(fastq_dir).parts
        sample_hint, read_key = _classify_fastq_name(fastq_file.name)

        if len(rel_parts) > 1:
            sample_acc = rel_parts[0]
        else:
            sample_acc = sample_hint or _strip_fastq_suffix(fastq_file.name)

        read_bucket = sample_to_reads.setdefault(
            sample_acc,
            {"R1": [], "R2": [], "single": []},
        )
        read_bucket[read_key].append(fastq_file)

    for read_bucket in sample_to_reads.values():
        for key in ("R1", "R2", "single"):
            read_bucket[key] = sorted(read_bucket[key])

    return sample_to_reads


def _merge_or_passthrough(files: List[Path], destination: Path) -> Path:
    """
    1 ファイルならそのまま返し、複数なら gzip バイト列を順に連結します。
    """
    if len(files) == 1:
        return files[0]

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as dst:
        for src in files:
            with src.open("rb") as fh:
                shutil.copyfileobj(fh, dst)
    return destination


def merge_lanes_by_cat(
    sample_to_reads: Dict[str, Dict[str, List[Path]]],
    merged_dir: Path,
    logger: logging.Logger,
) -> Dict[str, List[Path]]:
    """
    サンプルごとのサブレーン FASTQ をまとめ、解析用の FASTQ 一覧を返します。
    """
    merged_dir = Path(merged_dir)
    merged_dir.mkdir(parents=True, exist_ok=True)

    sample_to_fastqs: Dict[str, List[Path]] = {}

    for sample_acc, read_bucket in sorted(sample_to_reads.items()):
        r1_files = sorted(read_bucket.get("R1", []))
        r2_files = sorted(read_bucket.get("R2", []))
        single_files = sorted(read_bucket.get("single", []))

        if r1_files or r2_files:
            if single_files:
                logger.warning(
                    "paired-end と single-end FASTQ が混在しているため single を無視します: %s",
                    sample_acc,
                )

            if r1_files and r2_files:
                if len(r1_files) != len(r2_files):
                    logger.warning(
                        "R1/R2 の本数が一致しませんが、そのままマージします: %s (R1=%d, R2=%d)",
                        sample_acc,
                        len(r1_files),
                        len(r2_files),
                    )

                r1_merged = _merge_or_passthrough(
                    r1_files,
                    merged_dir / f"{sample_acc}.R1.fastq.gz",
                )
                r2_merged = _merge_or_passthrough(
                    r2_files,
                    merged_dir / f"{sample_acc}.R2.fastq.gz",
                )
                sample_to_fastqs[sample_acc] = [r1_merged, r2_merged]
                continue

            logger.warning(
                "片側だけの paired-end FASTQ を single-end として扱います: %s",
                sample_acc,
            )
            lone_reads = r1_files or r2_files
            sample_to_fastqs[sample_acc] = [
                _merge_or_passthrough(
                    lone_reads,
                    merged_dir / f"{sample_acc}.single.fastq.gz",
                )
            ]
            continue

        if single_files:
            sample_to_fastqs[sample_acc] = [
                _merge_or_passthrough(
                    single_files,
                    merged_dir / f"{sample_acc}.single.fastq.gz",
                )
            ]

    return sample_to_fastqs
