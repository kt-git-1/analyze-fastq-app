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
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class FastqRun:
    sample_acc: str
    run_id: str
    fastq_files: List[Path]
    rg_library: Optional[str] = None


@dataclass(frozen=True)
class _FastqNameParts:
    sample_hint: Optional[str]
    read_key: str
    run_hint: Optional[str]
    lane: Optional[str] = None


def _strip_fastq_suffix(filename: str) -> str:
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return filename


def _parse_fastq_name(filename: str) -> _FastqNameParts:
    stem = _strip_fastq_suffix(filename)

    # Illumina: SAMPLE_L001_R1_001
    m = re.match(
        r"^(?P<sample>.+?)_(?P<lane>L\d{3})_R(?P<read>[12])(?:_\d+)?$", stem
    )
    if m:
        sample = m.group("sample")
        lane = m.group("lane")
        return _FastqNameParts(sample, f"R{m.group('read')}", f"{sample}_{lane}", lane)

    # Filgen: Horse01_L1_1
    m = re.match(r"^(?P<sample>.+?)_(?P<lane>L\d+)_(?P<read>[12])$", stem)
    if m:
        sample = m.group("sample")
        lane = m.group("lane")
        return _FastqNameParts(sample, f"R{m.group('read')}", f"{sample}_{lane}", lane)

    # General R1/R2: sampleB_R1
    m = re.match(
        r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]?\d+)?$", stem, re.IGNORECASE
    )
    if m:
        return _FastqNameParts(m.group("sample"), f"R{m.group('read')}", None)

    # _1/_2: sampleC_1
    m = re.match(r"^(?P<sample>.+?)[._-](?P<read>[12])$", stem)
    if m:
        return _FastqNameParts(m.group("sample"), f"R{m.group('read')}", None)

    return _FastqNameParts(None, "single", None)


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
    parts = _parse_fastq_name(filename)
    return parts.sample_hint, parts.read_key, parts.run_hint


def _sanitize_run_id(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    return sanitized or "run"


def _infer_rg_library(filename: str) -> Optional[str]:
    stem = _strip_fastq_suffix(filename)
    match = re.search(r"(SCY\d+(?:\.\d+)?)", stem)
    if match:
        return match.group(1)
    return None


def _unique_run_id(sample_acc: str, base_run_id: str, used_run_ids: Dict[str, int]) -> str:
    base = _sanitize_run_id(base_run_id)
    key = f"{sample_acc}\0{base}"
    count = used_run_ids.get(key, 0) + 1
    used_run_ids[key] = count
    if count == 1:
        return base
    return f"{base}_{count}"


# ------------------------------------------------------------------
# ラン単位グルーピング（推奨）
# ------------------------------------------------------------------

def group_fastqs_by_run(
    fastq_dir: Path,
) -> Dict[str, List[FastqRun]]:
    """
    FASTQ ディレクトリからサンプル → ラン単位の FASTQ リストを返す。

    ENA 構成 (SAMPLE/RUN_DIR/file.gz) や Illumina レーン構成に対応する。
    各ランは独立にマッピングされ、あとでサンプル単位にマージされる前提。

    Returns
    -------
    dict[str, list[FastqRun]]
        {sample_acc: [FastqRun(...), ...]}
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
        parsed = _parse_fastq_name(fq.name)
        sample_hint = parsed.sample_hint
        read_key = parsed.read_key
        run_hint = parsed.run_hint

        is_run_folder_layout = (
            len(rel_parts) == 2
            and sample_hint is not None
            and parsed.lane is not None
            and rel_parts[0] != sample_hint
        )

        if is_run_folder_layout:
            sample_acc = sample_hint
        elif len(rel_parts) > 1:
            sample_acc = rel_parts[0]
        else:
            sample_acc = sample_hint or _strip_fastq_suffix(fq.name)

        if is_run_folder_layout:
            run_id = f"{rel_parts[0]}_{parsed.lane}"
        elif len(rel_parts) >= 3:
            run_id = rel_parts[-2]
        elif run_hint:
            run_id = run_hint
        else:
            run_id = sample_acc

        key: RunKey = (sample_acc, run_id)
        bucket = run_buckets.setdefault(key, {"R1": [], "R2": [], "single": []})
        bucket[read_key].append(fq)

    result: Dict[str, List[FastqRun]] = {}
    used_run_ids: Dict[str, int] = {}

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
            if len(r1) != len(r2):
                logger.warning(
                    "R1/R2 の本数が一致しないため余剰リードを除外: %s/%s (R1=%d, R2=%d)",
                    sample_acc,
                    run_id,
                    len(r1),
                    len(r2),
                )
            for idx, (r1_path, r2_path) in enumerate(zip(r1, r2), 1):
                base_run_id = run_id if idx == 1 else f"{run_id}_{idx}"
                unique_id = _unique_run_id(sample_acc, base_run_id, used_run_ids)
                rg_library = _infer_rg_library(r1_path.name) or _infer_rg_library(r2_path.name)
                result.setdefault(sample_acc, []).append(
                    FastqRun(sample_acc, unique_id, [r1_path, r2_path], rg_library)
                )
            continue

        if r1 or r2:
            lone_reads = r1 or r2
            logger.warning(
                "片側だけの paired-end FASTQ を除外します: %s/%s (%d ファイル)",
                sample_acc,
                run_id,
                len(lone_reads),
            )
            continue

        for single in singles:
            if run_id == sample_acc:
                base_run_id = _strip_fastq_suffix(single.name)
            else:
                base_run_id = run_id
            unique_id = _unique_run_id(sample_acc, base_run_id, used_run_ids)
            result.setdefault(sample_acc, []).append(
                FastqRun(
                    sample_acc,
                    unique_id,
                    [single],
                    _infer_rg_library(single.name),
                )
            )

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
