"""
FASTQ ファイルの検出・分類・ラン単位グルーピングを行うユーティリティ。

主なエントリポイント:
  - group_fastqs_by_run(fastq_dir) : ディレクトリ内の FASTQ をサンプル/ラン単位にグルーピング
  - parse_fastq_general(fastq_dir) : (後方互換) サンプル/R1/R2/single に分類
  - merge_lanes_by_cat(...)        : (後方互換) 複数レーンを1サンプル1セットに統合
"""

import logging
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def _strip_fastq_suffix(filename: str) -> str:
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def _classify_fastq_name(
    filename: str,
) -> Tuple[Optional[str], str, Optional[str]]:
    """
    FASTQ ファイル名から sample 名のヒント、read 種別、ラン（レーン）ヒントを推定する。

    Returns
    -------
    tuple[str | None, str, str | None]
        (sample_hint, read_key, run_hint)
    """
    stem = _strip_fastq_suffix(filename)

    # Illumina: SAMPLE_L001_R1_001
    m = re.match(
        r"^(?P<sample>.+?)_(?P<lane>L\d{3})_R(?P<read>[12])(?:_\d+)?$", stem
    )
    if m:
        return (
            m.group("sample"),
            f"R{m.group('read')}",
            f"{m.group('sample')}_{m.group('lane')}",
        )

    # Filgen: Horse01_L1_1
    m = re.match(r"^(?P<sample>.+?)_(?P<lane>L\d+)_(?P<read>[12])$", stem)
    if m:
        return (
            m.group("sample"),
            f"R{m.group('read')}",
            f"{m.group('sample')}_{m.group('lane')}",
        )

    # General R1/R2: sampleB_R1
    m = re.match(
        r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]?\d+)?$", stem, re.IGNORECASE
    )
    if m:
        return m.group("sample"), f"R{m.group('read')}", None

    # _1/_2: sampleC_1
    m = re.match(r"^(?P<sample>.+?)[._-](?P<read>[12])$", stem)
    if m:
        return m.group("sample"), f"R{m.group('read')}", None

    return None, "single", None


# ------------------------------------------------------------------
# ラン単位グルーピング（推奨）
# ------------------------------------------------------------------

def group_fastqs_by_run(
    fastq_dir: Path,
) -> Dict[str, List[Tuple[str, List[Path]]]]:
    """
    FASTQ ディレクトリからサンプル → ラン単位の FASTQ リストを返す。

    ENA 構成 (SAMPLE/RUN_DIR/file.gz) や Illumina レーン構成に対応する。
    各ランは独立にマッピングされ、あとでサンプル単位にマージされる前提。

    Returns
    -------
    dict[str, list[tuple[str, list[Path]]]]
        {sample_acc: [(run_id, [fastq_path, ...]), ...]}
        paired-end ラン: [R1_path, R2_path]
        single-end ラン: [single_path]
    """
    fastq_dir = Path(fastq_dir)

    RunKey = Tuple[str, str]
    run_buckets: Dict[RunKey, Dict[str, List[Path]]] = {}

    fastq_files = sorted(
        [
            *fastq_dir.rglob("*.fastq"),
            *fastq_dir.rglob("*.fastq.gz"),
            *fastq_dir.rglob("*.fq"),
            *fastq_dir.rglob("*.fq.gz"),
        ]
    )

    for fq in fastq_files:
        rel_parts = fq.relative_to(fastq_dir).parts
        sample_hint, read_key, run_hint = _classify_fastq_name(fq.name)

        if len(rel_parts) > 1:
            sample_acc = rel_parts[0]
        else:
            sample_acc = sample_hint or _strip_fastq_suffix(fq.name)

        if len(rel_parts) >= 3:
            run_id = rel_parts[-2]
        elif run_hint:
            run_id = run_hint
        else:
            run_id = sample_acc

        key: RunKey = (sample_acc, run_id)
        bucket = run_buckets.setdefault(key, {"R1": [], "R2": [], "single": []})
        bucket[read_key].append(fq)

    result: Dict[str, List[Tuple[str, List[Path]]]] = {}

    for (sample_acc, run_id), bucket in sorted(run_buckets.items()):
        r1 = sorted(bucket["R1"])
        r2 = sorted(bucket["R2"])
        singles = sorted(bucket["single"])

        if r1 and r2:
            if singles:
                logger.warning(
                    "paired-end と single-end が混在 (single を無視): %s/%s",
                    sample_acc,
                    run_id,
                )
            fastqs = [r1[0], r2[0]]
        elif r1:
            fastqs = r1[:1]
        elif r2:
            fastqs = r2[:1]
        elif singles:
            fastqs = singles[:1]
        else:
            continue

        result.setdefault(sample_acc, []).append((run_id, fastqs))

    return result


# ------------------------------------------------------------------
# 以下は後方互換のために残す（旧方式: 全ランを cat で結合）
# ------------------------------------------------------------------

def parse_fastq_general(fastq_dir: Path) -> Dict[str, Dict[str, List[Path]]]:
    """FASTQ を再帰探索して sample ごとに R1/R2/single に分類する（後方互換）。"""
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
        sample_hint, read_key, _ = _classify_fastq_name(fastq_file.name)

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
    log: logging.Logger,
) -> Dict[str, List[Path]]:
    """サンプルごとの FASTQ をまとめる（後方互換）。"""
    merged_dir = Path(merged_dir)
    merged_dir.mkdir(parents=True, exist_ok=True)
    sample_to_fastqs: Dict[str, List[Path]] = {}

    for sample_acc, read_bucket in sorted(sample_to_reads.items()):
        r1_files = sorted(read_bucket.get("R1", []))
        r2_files = sorted(read_bucket.get("R2", []))
        single_files = sorted(read_bucket.get("single", []))

        if r1_files or r2_files:
            if single_files:
                log.warning(
                    "paired-end と single-end FASTQ が混在しているため single を無視します: %s",
                    sample_acc,
                )
            if r1_files and r2_files:
                if len(r1_files) != len(r2_files):
                    log.warning(
                        "R1/R2 の本数が一致しませんが、そのままマージします: %s (R1=%d, R2=%d)",
                        sample_acc, len(r1_files), len(r2_files),
                    )
                sample_to_fastqs[sample_acc] = [
                    _merge_or_passthrough(r1_files, merged_dir / f"{sample_acc}.R1.fastq.gz"),
                    _merge_or_passthrough(r2_files, merged_dir / f"{sample_acc}.R2.fastq.gz"),
                ]
                continue
            log.warning("片側だけの paired-end FASTQ を single-end として扱います: %s", sample_acc)
            lone_reads = r1_files or r2_files
            sample_to_fastqs[sample_acc] = [
                _merge_or_passthrough(lone_reads, merged_dir / f"{sample_acc}.single.fastq.gz")
            ]
            continue

        if single_files:
            sample_to_fastqs[sample_acc] = [
                _merge_or_passthrough(single_files, merged_dir / f"{sample_acc}.single.fastq.gz")
            ]

    return sample_to_fastqs
