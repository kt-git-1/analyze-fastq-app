import csv
import logging
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PCASite:
    chrom: str
    pos: int
    site_id: str
    ref: str = ""
    alt: str = ""


@dataclass(frozen=True)
class PCAQCConfig:
    min_mapq: int = 30
    min_baseq: int = 30
    trim_ends: int = 2
    max_sample_missing: float = 0.9
    max_site_missing: float = 0.9
    min_maf: float = 0.0
    exclude_sex_chr: bool = False
    ld_window: int = 50
    ld_step: int = 5
    ld_r2: float = 0.2


@dataclass(frozen=True)
class MatrixStats:
    input_samples: int
    input_sites: int
    kept_samples: int
    kept_sites: int
    missingness_removed_sites: int
    monomorphic_removed_sites: int
    maf_removed_sites: int
    sex_chr_removed_sites: int


def _is_nonempty_file(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _stage_done(path: Path, force: bool) -> bool:
    return (not force) and _is_nonempty_file(path)


def _stages_done(paths: Iterable[Path], force: bool) -> bool:
    return (not force) and all(_is_nonempty_file(path) for path in paths)


def _run_command(cmd: List[str], *, cwd: Optional[Path] = None) -> None:
    logger.info("外部コマンドを実行します: %s", " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def _final_dedup_bam(config, sample: str) -> Path:
    return config.results_dir / sample / "dedup" / ("%s.dedup.sorted.bam" % sample)


def _pca_qc_from_args(args) -> PCAQCConfig:
    return PCAQCConfig(
        min_mapq=getattr(args, "pca_min_mapq", 30),
        min_baseq=getattr(args, "pca_min_baseq", 30),
        trim_ends=getattr(args, "pca_trim_ends", 2),
        max_sample_missing=getattr(args, "pca_max_sample_missing", 0.9),
        max_site_missing=getattr(args, "pca_max_site_missing", 0.9),
        min_maf=getattr(args, "pca_min_maf", 0.0),
        exclude_sex_chr=getattr(args, "pca_exclude_sex_chr", False),
        ld_window=getattr(args, "pca_ld_window", 50),
        ld_step=getattr(args, "pca_ld_step", 5),
        ld_r2=getattr(args, "pca_ld_r2", 0.2),
    )


def _median_mapped_read_length(bam_path: Path, max_reads: int = 20000) -> Optional[float]:
    if not _is_nonempty_file(bam_path):
        return None
    import pysam

    lengths: List[int] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.query_length is None:
                continue
            lengths.append(read.query_length)
            if len(lengths) >= max_reads:
                break
    if not lengths:
        return None
    return float(statistics.median(lengths))


def _read_auto_inferred_type(summary: Path) -> Optional[str]:
    if not _is_nonempty_file(summary):
        return None
    with summary.open(errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) >= 2 and row[0] == "inferred_data_type" and row[1] in {"ancient", "modern"}:
                return row[1]
    return None


def infer_pca_data_type(config, samples: Iterable[str], cohort_dir: Path, *, force: bool = False) -> str:
    requested = config.data_type
    if requested in {"ancient", "modern"}:
        return requested
    if requested != "auto":
        raise NotImplementedError("PCA stage supports data_type=ancient, modern, or auto")

    summary = cohort_dir / "auto_data_type_summary.tsv"
    existing = None if force else _read_auto_inferred_type(summary)
    if existing:
        logger.info("再開: 既存auto data_type推定を使用します: %s (%s)", existing, summary)
        return existing

    rows: List[Tuple[str, str, str]] = []
    medians: List[float] = []
    for sample in sorted(set(samples)):
        bam_path = _final_dedup_bam(config, sample)
        try:
            median_len = _median_mapped_read_length(bam_path)
        except Exception as exc:
            logger.warning("auto data_type 推定でBAM read lengthを読めませんでした: %s (%s)", sample, exc)
            median_len = None
        if median_len is None:
            rows.append([sample, "", "no_mapped_read_length"])
        else:
            medians.append(median_len)
            rows.append([sample, "%.3f" % median_len, "ok"])

    if not medians:
        raise ValueError("data_type=auto ですが、final dedup BAMからread lengthを推定できません")

    cohort_median = float(statistics.median(medians))
    inferred = "ancient" if cohort_median <= 100.0 else "modern"
    summary.parent.mkdir(parents=True, exist_ok=True)
    with summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["requested_data_type", requested])
        writer.writerow(["inferred_data_type", inferred])
        writer.writerow(["cohort_median_read_length", "%.3f" % cohort_median])
        writer.writerow(["read_length_cutoff", "100.000"])
        writer.writerow([])
        writer.writerow(["sample", "median_mapped_read_length", "status"])
        writer.writerows(rows)
    logger.info(
        "data_type=auto: cohort median read length %.3f bp から %s と推定しました: %s",
        cohort_median,
        inferred,
        summary,
    )
    return inferred


def read_pca_sites(path: Path) -> List[PCASite]:
    """Read target SNP sites from VCF or BED-like text."""
    if not path.exists():
        raise FileNotFoundError("PCA sites file not found: %s" % path)

    sites: List[PCASite] = []
    with path.open(errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if not line.strip() or line.startswith("##"):
                continue
            if line.startswith("site_id\t"):
                continue
            if line.startswith("#CHROM"):
                continue

            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) >= 3 and path.name == "pca_sites.tsv":
                site_id = parts[0]
                chrom = parts[1]
                pos_text = parts[2]
                ref = parts[3] if len(parts) >= 4 else ""
                alt = parts[4] if len(parts) >= 5 else ""
                sites.append(PCASite(chrom, int(pos_text), site_id, ref.upper(), alt.upper()))
                continue

            if len(parts) >= 5 and not path.suffix.lower().startswith(".bed"):
                chrom, pos_text, site_id, ref, alt = parts[:5]
                if "," in alt or len(ref) != 1 or len(alt) != 1:
                    continue
                pos = int(pos_text)
                if site_id == ".":
                    site_id = "%s:%d:%s:%s" % (chrom, pos, ref, alt)
                sites.append(PCASite(chrom, pos, site_id, ref.upper(), alt.upper()))
                continue

            if len(parts) < 3:
                continue
            chrom = parts[0]
            # BED uses 0-based start. Store 1-based genomic position.
            pos = int(parts[1]) + 1
            site_id = parts[3] if len(parts) >= 4 else "%s:%d" % (chrom, pos)
            ref = parts[4].upper() if len(parts) >= 5 else ""
            alt = parts[5].upper() if len(parts) >= 6 else ""
            if ref and alt and (len(ref) != 1 or len(alt) != 1):
                continue
            sites.append(PCASite(chrom, pos, site_id, ref, alt))

    if not sites:
        raise ValueError("PCA sites file contains no usable biallelic SNP sites: %s" % path)
    return sites


def write_sites_table(sites: List[PCASite], out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["site_id", "chrom", "pos", "ref", "alt"])
        for site in sites:
            writer.writerow([site.site_id, site.chrom, site.pos, site.ref, site.alt])
    return out_path


def _choose_pileup_allele(
    column,
    site: PCASite,
    min_baseq: int,
    trim_ends: int,
) -> Tuple[str, int, int]:
    candidates: List[Tuple[int, int, str]] = []
    allowed = {base for base in (site.ref, site.alt) if base}
    if not allowed:
        allowed = {"A", "C", "G", "T"}

    for pileup_read in column.pileups:
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        read = pileup_read.alignment
        qpos = pileup_read.query_position
        if qpos is None or read.query_sequence is None:
            continue
        if trim_ends and (qpos < trim_ends or qpos >= read.query_length - trim_ends):
            continue
        base = read.query_sequence[qpos].upper()
        if base not in allowed:
            continue
        baseq = read.query_qualities[qpos] if read.query_qualities is not None else 0
        if baseq < min_baseq:
            continue
        candidates.append((baseq, read.mapping_quality, base))

    if not candidates:
        return "", 0, 0
    candidates.sort(reverse=True)
    baseq, mapq, base = candidates[0]
    return base, baseq, mapq


def extract_pseudohaploid_calls(
    sample_bams: Dict[str, Path],
    sites: List[PCASite],
    out_path: Path,
    *,
    min_mapq: int = 30,
    min_baseq: int = 30,
    trim_ends: int = 2,
) -> Path:
    """Extract one high-quality allele per sample/site from final dedup BAMs."""
    import pysam

    sites_by_chrom: Dict[str, List[Tuple[int, PCASite]]] = {}
    for idx, site in enumerate(sites):
        sites_by_chrom.setdefault(site.chrom, []).append((idx, site))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "site_id", "chrom", "pos", "ref", "alt", "allele", "baseq", "mapq"])
        for sample, bam_path in sorted(sample_bams.items()):
            if not _is_nonempty_file(bam_path):
                logger.warning("PCA allele extraction skipped missing BAM: %s (%s)", sample, bam_path)
                continue
            logger.info("pseudo-haploid allele を抽出します: %s", sample)
            calls_by_site: Dict[int, Tuple[str, int, int]] = {}
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                references = set(bam.references)
                for chrom, chrom_sites in sites_by_chrom.items():
                    if chrom not in references:
                        logger.warning("PCA sites contig is absent from BAM header: %s (%s)", chrom, sample)
                        continue
                    positions = {}
                    for idx, site in chrom_sites:
                        positions.setdefault(site.pos - 1, []).append((idx, site))
                    start = min(positions)
                    stop = max(positions) + 1
                    for column in bam.pileup(
                        chrom,
                        start,
                        stop,
                        truncate=True,
                        stepper="samtools",
                        min_mapping_quality=min_mapq,
                    ):
                        target_sites = positions.get(column.reference_pos)
                        if not target_sites:
                            continue
                        for idx, site in target_sites:
                            calls_by_site[idx] = _choose_pileup_allele(column, site, min_baseq, trim_ends)
                for idx, site in enumerate(sites):
                    allele, baseq, mapq = calls_by_site.get(idx, ("", 0, 0))
                    writer.writerow([sample, site.site_id, site.chrom, site.pos, site.ref, site.alt, allele, baseq, mapq])
    return out_path


def _load_raw_calls(path: Path) -> Dict[str, Dict[str, str]]:
    calls: Dict[str, Dict[str, str]] = {}
    with path.open(errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            site_id = row["site_id"]
            calls.setdefault(sample, {})[site_id] = row.get("allele", "")
    return calls


def build_pseudohaploid_matrix(
    raw_calls_path: Path,
    sites: List[PCASite],
    samples: List[str],
    out_path: Path,
) -> Path:
    calls = _load_raw_calls(raw_calls_path)
    site_lookup = {site.site_id: site for site in sites}
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + [site.site_id for site in sites])
        for sample in sorted(samples):
            row = [sample]
            for site in sites:
                allele = calls.get(sample, {}).get(site.site_id, "")
                code = ""
                if allele:
                    if site.ref and site.alt:
                        if allele == site.ref:
                            code = "0"
                        elif allele == site.alt:
                            code = "1"
                    else:
                        observed = [
                            calls.get(s, {}).get(site.site_id, "")
                            for s in samples
                            if calls.get(s, {}).get(site.site_id, "")
                        ]
                        if observed:
                            major = sorted(set(observed), key=lambda b: (-observed.count(b), b))[0]
                            code = "0" if allele == major else "1"
                row.append(code)
            writer.writerow(row)
    # Validate every requested site existed in the lookup; keeps linters quiet and catches accidental duplicates.
    if len(site_lookup) != len(sites):
        logger.warning("PCA sites include duplicate IDs; downstream matrix columns may be ambiguous")
    return out_path


def _count_pileup_alleles(
    column,
    site: PCASite,
    min_baseq: int,
    trim_ends: int,
) -> Tuple[int, int, int]:
    ref_count = 0
    alt_count = 0
    depth = 0
    for pileup_read in column.pileups:
        if pileup_read.is_del or pileup_read.is_refskip:
            continue
        read = pileup_read.alignment
        qpos = pileup_read.query_position
        if qpos is None or read.query_sequence is None:
            continue
        if trim_ends and (qpos < trim_ends or qpos >= read.query_length - trim_ends):
            continue
        baseq = read.query_qualities[qpos] if read.query_qualities is not None else 0
        if baseq < min_baseq:
            continue
        base = read.query_sequence[qpos].upper()
        if site.ref and base == site.ref:
            ref_count += 1
            depth += 1
        elif site.alt and base == site.alt:
            alt_count += 1
            depth += 1
    return ref_count, alt_count, depth


def _diploid_dosage_from_counts(ref_count: int, alt_count: int, depth: int) -> str:
    if depth <= 0:
        return ""
    alt_fraction = alt_count / depth
    if alt_fraction <= 0.2:
        return "0"
    if alt_fraction >= 0.8:
        return "2"
    return "1"


def extract_modern_diploid_calls(
    sample_bams: Dict[str, Path],
    sites: List[PCASite],
    out_path: Path,
    *,
    min_mapq: int = 30,
    min_baseq: int = 30,
    trim_ends: int = 2,
) -> Path:
    """Extract diploid ref/alt dosage per sample/site from final dedup BAMs."""
    import pysam

    sites_by_chrom: Dict[str, List[Tuple[int, PCASite]]] = {}
    for idx, site in enumerate(sites):
        sites_by_chrom.setdefault(site.chrom, []).append((idx, site))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "site_id", "chrom", "pos", "ref", "alt", "ref_count", "alt_count", "depth", "dosage"])
        for sample, bam_path in sorted(sample_bams.items()):
            calls_by_site: Dict[int, Tuple[int, int, int, str]] = {}
            if not _is_nonempty_file(bam_path):
                logger.warning("modern PCA genotype extraction skipped missing BAM: %s (%s)", sample, bam_path)
            else:
                logger.info("modern diploid genotype をBAMから抽出します: %s", sample)
                with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                    references = set(bam.references)
                    for chrom, chrom_sites in sites_by_chrom.items():
                        if chrom not in references:
                            logger.warning("PCA sites contig is absent from BAM header: %s (%s)", chrom, sample)
                            continue
                        positions = {}
                        for idx, site in chrom_sites:
                            positions.setdefault(site.pos - 1, []).append((idx, site))
                        start = min(positions)
                        stop = max(positions) + 1
                        for column in bam.pileup(
                            chrom,
                            start,
                            stop,
                            truncate=True,
                            stepper="samtools",
                            min_mapping_quality=min_mapq,
                        ):
                            target_sites = positions.get(column.reference_pos)
                            if not target_sites:
                                continue
                            for idx, site in target_sites:
                                ref_count, alt_count, depth = _count_pileup_alleles(column, site, min_baseq, trim_ends)
                                dosage = _diploid_dosage_from_counts(ref_count, alt_count, depth)
                                calls_by_site[idx] = (ref_count, alt_count, depth, dosage)
            for idx, site in enumerate(sites):
                ref_count, alt_count, depth, dosage = calls_by_site.get(idx, (0, 0, 0, ""))
                writer.writerow([sample, site.site_id, site.chrom, site.pos, site.ref, site.alt, ref_count, alt_count, depth, dosage])
    return out_path


def build_modern_genotype_matrix(
    genotype_calls_path: Path,
    sites: List[PCASite],
    samples: List[str],
    out_path: Path,
) -> Path:
    calls: Dict[str, Dict[str, str]] = {}
    with genotype_calls_path.open(errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            calls.setdefault(row["sample"], {})[row["site_id"]] = row.get("dosage", "")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + [site.site_id for site in sites])
        for sample in sorted(samples):
            writer.writerow([sample] + [calls.get(sample, {}).get(site.site_id, "") for site in sites])
    return out_path


def filter_matrix(
    matrix_path: Path,
    out_path: Path,
    *,
    max_site_missing: float = 0.9,
    max_sample_missing: float = 0.9,
    min_maf: float = 0.0,
    exclude_sex_chr: bool = False,
    sites: Optional[List[PCASite]] = None,
    ploidy: int = 1,
) -> Path:
    _filter_matrix(matrix_path, out_path, max_site_missing, max_sample_missing, min_maf, exclude_sex_chr, sites, ploidy)
    return out_path


def _is_sex_chrom(chrom: str) -> bool:
    normalized = chrom.lower().replace("chr", "")
    return normalized in {"x", "y", "w", "z", "23", "24"}


def _filter_matrix(
    matrix_path: Path,
    out_path: Path,
    max_site_missing: float,
    max_sample_missing: float,
    min_maf: float,
    exclude_sex_chr: bool,
    sites: Optional[List[PCASite]],
    ploidy: int,
) -> MatrixStats:
    matrix_stats, kept_rows, kept_site_ids, data = _filter_matrix_components(
        matrix_path,
        max_site_missing,
        max_sample_missing,
        min_maf,
        exclude_sex_chr,
        sites,
        ploidy,
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + kept_site_ids)
        for row, values in zip(kept_rows, data):
            writer.writerow([row[0]] + ["" if np.isnan(value) else str(int(value)) for value in values])

    return matrix_stats


def _filter_matrix_stats(
    matrix_path: Path,
    max_site_missing: float,
    max_sample_missing: float,
    min_maf: float,
    exclude_sex_chr: bool,
    sites: Optional[List[PCASite]],
    ploidy: int,
) -> MatrixStats:
    matrix_stats, _kept_rows, _kept_site_ids, _data = _filter_matrix_components(
        matrix_path,
        max_site_missing,
        max_sample_missing,
        min_maf,
        exclude_sex_chr,
        sites,
        ploidy,
    )
    return matrix_stats


def _filter_matrix_components(
    matrix_path: Path,
    max_site_missing: float,
    max_sample_missing: float,
    min_maf: float,
    exclude_sex_chr: bool,
    sites: Optional[List[PCASite]],
    ploidy: int,
) -> Tuple[MatrixStats, List[List[str]], List[str], np.ndarray]:
    with matrix_path.open(errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)

    if not rows or len(header) < 2:
        raise ValueError("PCA matrix has no sample or site data: %s" % matrix_path)

    data = np.array([[np.nan if value == "" else float(value) for value in row[1:]] for row in rows], dtype=float)
    sample_missing = np.mean(np.isnan(data), axis=1)
    keep_samples = sample_missing <= max_sample_missing
    data = data[keep_samples, :]
    kept_rows = [row for row, keep in zip(rows, keep_samples) if keep]

    if data.size == 0 or not kept_rows:
        raise ValueError("No samples remain after PCA missingness filtering")

    site_missing = np.mean(np.isnan(data), axis=0)
    variable: List[bool] = []
    maf_keep: List[bool] = []
    for col_idx in range(data.shape[1]):
        values = data[:, col_idx]
        observed = values[~np.isnan(values)]
        variable.append(len(set(observed.tolist())) > 1)
        if len(observed) == 0:
            maf_keep.append(False)
        else:
            allele_freq = float(np.mean(observed) / max(1, ploidy))
            maf_keep.append(min(allele_freq, 1.0 - allele_freq) >= min_maf)

    site_lookup = {site.site_id: site for site in sites or []}
    sex_chr_keep: List[bool] = []
    for site_id in header[1:]:
        site = site_lookup.get(site_id)
        sex_chr_keep.append(not (exclude_sex_chr and site and _is_sex_chrom(site.chrom)))

    missing_keep = site_missing <= max_site_missing
    variable_keep = np.array(variable, dtype=bool)
    maf_keep_array = np.array(maf_keep, dtype=bool)
    sex_chr_keep_array = np.array(sex_chr_keep, dtype=bool)
    keep_sites = missing_keep & variable_keep & maf_keep_array & sex_chr_keep_array
    data = data[:, keep_sites]
    kept_site_ids = [site_id for site_id, keep in zip(header[1:], keep_sites) if keep]

    if data.size == 0 or not kept_site_ids:
        raise ValueError("No variable SNP sites remain after PCA filtering")

    matrix_stats = MatrixStats(
        input_samples=len(rows),
        input_sites=len(header) - 1,
        kept_samples=len(kept_rows),
        kept_sites=len(kept_site_ids),
        missingness_removed_sites=int(np.sum(~missing_keep)),
        monomorphic_removed_sites=int(np.sum(missing_keep & ~variable_keep)),
        maf_removed_sites=int(np.sum(missing_keep & variable_keep & ~maf_keep_array)),
        sex_chr_removed_sites=int(np.sum(missing_keep & variable_keep & maf_keep_array & ~sex_chr_keep_array)),
    )
    return matrix_stats, kept_rows, kept_site_ids, data


def _load_numeric_matrix(matrix_path: Path) -> Tuple[List[str], List[str], np.ndarray]:
    with matrix_path.open(errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    samples = [row[0] for row in rows]
    site_ids = header[1:]
    data = np.array([[np.nan if value == "" else float(value) for value in row[1:]] for row in rows], dtype=float)
    if len(samples) < 2:
        raise ValueError("PCA requires at least two samples")
    if data.shape[1] < 1:
        raise ValueError("PCA requires at least one variable SNP")
    return samples, site_ids, data


def _impute_and_scale(data: np.ndarray) -> np.ndarray:
    col_means = np.nanmean(data, axis=0)
    filled = np.where(np.isnan(data), col_means, data)
    centered = filled - np.mean(filled, axis=0)
    std = np.std(centered, axis=0)
    std[std == 0] = 1.0
    return centered / std


def run_pca_and_mds(matrix_path: Path, pca_dir: Path) -> Tuple[Path, Path, Path]:
    samples, _site_ids, data = _load_numeric_matrix(matrix_path)
    scaled = _impute_and_scale(data)
    u, singular_values, _vt = np.linalg.svd(scaled, full_matrices=False)
    components = min(10, u.shape[1])
    scores = u[:, :components] * singular_values[:components]
    eigenvalues = (singular_values ** 2) / max(1, scaled.shape[0] - 1)
    explained = eigenvalues / np.sum(eigenvalues) if np.sum(eigenvalues) else eigenvalues

    pca_dir.mkdir(parents=True, exist_ok=True)
    scores_path = pca_dir / "pca_scores.tsv"
    variance_path = pca_dir / "pca_variance.tsv"
    with scores_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + ["PC%d" % (idx + 1) for idx in range(components)])
        for sample, values in zip(samples, scores[:, :components]):
            writer.writerow([sample] + ["%.8g" % value for value in values])
    with variance_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["component", "eigenvalue", "explained_variance"])
        for idx, value in enumerate(eigenvalues[:components]):
            writer.writerow(["PC%d" % (idx + 1), "%.8g" % value, "%.8g" % explained[idx]])

    mds_path = pca_dir / "mds.tsv"
    distances = _pairwise_distances(scaled)
    mds_coords = _classical_mds(distances, dimensions=2)
    with mds_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "MDS1", "MDS2"])
        for sample, values in zip(samples, mds_coords):
            writer.writerow([sample, "%.8g" % values[0], "%.8g" % values[1]])
    return scores_path, variance_path, mds_path


def write_mds_from_matrix(matrix_path: Path, out_path: Path) -> Path:
    samples, _site_ids, data = _load_numeric_matrix(matrix_path)
    scaled = _impute_and_scale(data)
    distances = _pairwise_distances(scaled)
    mds_coords = _classical_mds(distances, dimensions=2)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "MDS1", "MDS2"])
        for sample, values in zip(samples, mds_coords):
            writer.writerow([sample, "%.8g" % values[0], "%.8g" % values[1]])
    return out_path


def _pairwise_distances(data: np.ndarray) -> np.ndarray:
    n_samples = data.shape[0]
    distances = np.zeros((n_samples, n_samples), dtype=float)
    for i in range(n_samples):
        diff = data[i] - data
        distances[i, :] = np.sqrt(np.sum(diff * diff, axis=1))
    return distances


def _classical_mds(distances: np.ndarray, dimensions: int = 2) -> np.ndarray:
    n_samples = distances.shape[0]
    if n_samples == 0:
        return np.empty((0, dimensions))
    squared = distances ** 2
    centering = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    gram = -0.5 * centering.dot(squared).dot(centering)
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    coords = np.zeros((n_samples, dimensions), dtype=float)
    usable = min(dimensions, len(eigenvalues))
    positive = np.maximum(eigenvalues[:usable], 0)
    coords[:, :usable] = eigenvectors[:, :usable] * np.sqrt(positive)
    return coords


def export_eigenstrat(matrix_path: Path, sites: List[PCASite], out_dir: Path) -> Tuple[Path, Path, Path]:
    samples, site_ids, data = _load_numeric_matrix(matrix_path)
    site_lookup = {site.site_id: site for site in sites}
    out_dir.mkdir(parents=True, exist_ok=True)
    geno_path = out_dir / "cohort.geno"
    snp_path = out_dir / "cohort.snp"
    ind_path = out_dir / "cohort.ind"

    with geno_path.open("w") as handle:
        for col_idx in range(data.shape[1]):
            values = []
            for value in data[:, col_idx]:
                values.append("9" if np.isnan(value) else str(int(value)))
            handle.write("".join(values) + "\n")
    with snp_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for site_id in site_ids:
            site = site_lookup.get(site_id, PCASite("NA", 0, site_id))
            writer.writerow([site.site_id, site.chrom, "0.0", site.pos, site.ref or "0", site.alt or "1"])
    with ind_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for sample in samples:
            writer.writerow([sample, "U", "Unknown"])
    return geno_path, snp_path, ind_path


def export_plink_text(
    matrix_path: Path,
    sites: List[PCASite],
    out_dir: Path,
    *,
    pseudo_haploid: bool = False,
) -> Tuple[Path, Path]:
    samples, site_ids, data = _load_numeric_matrix(matrix_path)
    site_lookup = {site.site_id: site for site in sites}
    out_dir.mkdir(parents=True, exist_ok=True)
    tped_path = out_dir / "cohort.tped"
    tfam_path = out_dir / "cohort.tfam"

    with tped_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for col_idx, site_id in enumerate(site_ids):
            site = site_lookup.get(site_id, PCASite("0", 0, site_id, "0", "1"))
            row = [site.chrom.replace("chr", ""), site.site_id, "0", str(site.pos)]
            for value in data[:, col_idx]:
                if np.isnan(value):
                    row.extend(["0", "0"])
                elif int(value) == 0:
                    row.extend([site.ref or "A", site.ref or "A"])
                elif pseudo_haploid:
                    row.extend([site.alt or "B", site.alt or "B"])
                elif int(value) == 1:
                    row.extend([site.ref or "A", site.alt or "B"])
                else:
                    row.extend([site.alt or "B", site.alt or "B"])
            writer.writerow(row)
    with tfam_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=" ")
        for sample in samples:
            writer.writerow([sample, sample, "0", "0", "0", "-9"])
    return tped_path, tfam_path


def _write_parfile(path: Path, values: Dict[str, Union[Path, str, int, float]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for key, value in values.items():
            handle.write("%s: %s\n" % (key, value))
    return path


def _convert_evec_to_tsv(evec_path: Path, out_path: Path) -> Path:
    rows: List[List[str]] = []
    components = 0
    with evec_path.open(errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            sample = parts[0]
            pcs = parts[1:-1]
            components = max(components, len(pcs))
            rows.append([sample] + pcs)

    if not rows:
        raise ValueError("smartpca evec file contains no sample scores: %s" % evec_path)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + ["PC%d" % (idx + 1) for idx in range(components)])
        for row in rows:
            writer.writerow(row)
    return out_path


def _convert_eval_to_tsv(eval_path: Path, out_path: Path) -> Path:
    eigenvalues: List[float] = []
    with eval_path.open(errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            try:
                eigenvalues.append(float(stripped.split()[0]))
            except ValueError:
                continue

    if not eigenvalues:
        raise ValueError("smartpca eval file contains no eigenvalues: %s" % eval_path)

    total = sum(eigenvalues)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["component", "eigenvalue", "explained_variance"])
        for idx, value in enumerate(eigenvalues):
            explained = value / total if total else 0.0
            writer.writerow(["PC%d" % (idx + 1), "%.8g" % value, "%.8g" % explained])
    return out_path


def _count_pruned_snps(prune_in: Path) -> int:
    if not prune_in.exists():
        return 0
    with prune_in.open(errors="replace") as handle:
        return sum(1 for line in handle if line.strip())


def run_eigensoft_pca(
    filtered_matrix: Path,
    sites: List[PCASite],
    cohort_dir: Path,
    qc: PCAQCConfig,
    *,
    pseudo_haploid: bool = False,
    force: bool = False,
) -> Tuple[Tuple[Path, Path, Path], Tuple[Path, Path], Path, Path, Path, int]:
    plink_dir = cohort_dir / "plink"
    eigenstrat_dir = cohort_dir / "eigenstrat"
    pca_dir = cohort_dir / "pca"
    par_dir = cohort_dir / "eigensoft"

    tped_path = plink_dir / "cohort.tped"
    tfam_path = plink_dir / "cohort.tfam"
    if _stages_done([tped_path, tfam_path], force):
        logger.info("再開: 既存PLINK text出力を使用します: %s", plink_dir)
        plink_files = (tped_path, tfam_path)
    else:
        plink_files = export_plink_text(filtered_matrix, sites, plink_dir, pseudo_haploid=pseudo_haploid)
    tfile_prefix = plink_dir / "cohort"
    bed_prefix = plink_dir / "cohort"
    qc_prefix = plink_dir / "cohort.qc"
    prune_prefix = plink_dir / "cohort.prune"
    pruned_prefix = plink_dir / "cohort.pruned"

    bed_files = [Path(str(bed_prefix) + suffix) for suffix in (".bed", ".bim", ".fam")]
    if _stages_done(bed_files, force):
        logger.info("再開: 既存PLINK binary出力を使用します: %s", bed_prefix)
    else:
        _run_command(["plink", "--tfile", str(tfile_prefix), "--make-bed", "--out", str(bed_prefix)])

    qc_files = [Path(str(qc_prefix) + suffix) for suffix in (".bed", ".bim", ".fam")]
    if _stages_done(qc_files, force):
        logger.info("再開: 既存PLINK QC出力を使用します: %s", qc_prefix)
    else:
        _run_command(
            [
                "plink",
                "--bfile",
                str(bed_prefix),
                "--mind",
                str(qc.max_sample_missing),
                "--geno",
                str(qc.max_site_missing),
                "--maf",
                str(qc.min_maf),
                "--make-bed",
                "--out",
                str(qc_prefix),
            ]
        )

    prune_in = Path(str(prune_prefix) + ".prune.in")
    if _stage_done(prune_in, force):
        logger.info("再開: 既存PLINK LD pruningリストを使用します: %s", prune_in)
    else:
        _run_command(
            [
                "plink",
                "--bfile",
                str(qc_prefix),
                "--indep-pairwise",
                str(qc.ld_window),
                str(qc.ld_step),
                str(qc.ld_r2),
                "--out",
                str(prune_prefix),
            ]
        )

    pruned_files = [Path(str(pruned_prefix) + suffix) for suffix in (".bed", ".bim", ".fam")]
    if _stages_done(pruned_files, force):
        logger.info("再開: 既存PLINK LD pruned出力を使用します: %s", pruned_prefix)
    else:
        _run_command(
            [
                "plink",
                "--bfile",
                str(qc_prefix),
                "--extract",
                str(prune_in),
                "--make-bed",
                "--out",
                str(pruned_prefix),
            ]
        )

    geno_path = eigenstrat_dir / "cohort.geno"
    snp_path = eigenstrat_dir / "cohort.snp"
    ind_path = eigenstrat_dir / "cohort.ind"
    eigenstrat_dir.mkdir(parents=True, exist_ok=True)
    convertf_par = _write_parfile(
        par_dir / "convertf.par",
        {
            "inputformat": "PACKEDPED",
            "genotypename": str(pruned_prefix) + ".bed",
            "snpname": str(pruned_prefix) + ".bim",
            "indivname": str(pruned_prefix) + ".fam",
            "outputformat": "EIGENSTRAT",
            "genotypeoutname": geno_path,
            "snpoutname": snp_path,
            "indivoutname": ind_path,
        },
    )
    if _stages_done([geno_path, snp_path, ind_path], force):
        logger.info("再開: 既存EIGENSTRAT出力を使用します: %s", eigenstrat_dir)
    else:
        _run_command(["convertf", "-p", str(convertf_par)])

    pca_dir.mkdir(parents=True, exist_ok=True)
    evec_path = pca_dir / "cohort.evec"
    eval_path = pca_dir / "cohort.eval"
    smartpca_par = _write_parfile(
        par_dir / "smartpca.par",
        {
            "genotypename": geno_path,
            "snpname": snp_path,
            "indivname": ind_path,
            "evecoutname": evec_path,
            "evaloutname": eval_path,
            "numoutevec": 10,
        },
    )
    if _stages_done([evec_path, eval_path], force):
        logger.info("再開: 既存smartpca出力を使用します: %s", pca_dir)
    else:
        _run_command(["smartpca", "-p", str(smartpca_par)])

    pca_scores = pca_dir / "pca_scores.tsv"
    if _stage_done(pca_scores, force):
        logger.info("再開: 既存PCA score TSVを使用します: %s", pca_scores)
    else:
        pca_scores = _convert_evec_to_tsv(evec_path, pca_scores)
    pca_variance = pca_dir / "pca_variance.tsv"
    if _stage_done(pca_variance, force):
        logger.info("再開: 既存PCA variance TSVを使用します: %s", pca_variance)
    else:
        pca_variance = _convert_eval_to_tsv(eval_path, pca_variance)
    mds_path = pca_dir / "mds.tsv"
    if _stage_done(mds_path, force):
        logger.info("再開: 既存MDS TSVを使用します: %s", mds_path)
    else:
        mds_path = write_mds_from_matrix(filtered_matrix, mds_path)

    return (
        (geno_path, snp_path, ind_path),
        plink_files,
        pca_scores,
        pca_variance,
        mds_path,
        _count_pruned_snps(prune_in),
    )


def _write_qc_summary(
    out_path: Path,
    samples: List[str],
    sites: List[PCASite],
    filtered_matrix: Path,
    qc: PCAQCConfig,
    engine: str,
    matrix_stats: MatrixStats,
    data_type: str,
    ploidy: int,
    requested_data_type: Optional[str] = None,
    ld_pruned_sites: Optional[int] = None,
) -> Path:
    kept_samples, kept_sites, _data = _load_numeric_matrix(filtered_matrix)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        if requested_data_type:
            writer.writerow(["requested_data_type", requested_data_type])
        writer.writerow(["data_type", data_type])
        writer.writerow(["ploidy", ploidy])
        writer.writerow(["engine", engine])
        writer.writerow(["min_mapq", qc.min_mapq])
        writer.writerow(["min_baseq", qc.min_baseq])
        writer.writerow(["trim_ends", qc.trim_ends])
        writer.writerow(["max_sample_missing", qc.max_sample_missing])
        writer.writerow(["max_site_missing", qc.max_site_missing])
        writer.writerow(["min_maf", qc.min_maf])
        writer.writerow(["exclude_sex_chr", qc.exclude_sex_chr])
        writer.writerow(["ld_window", qc.ld_window])
        writer.writerow(["ld_step", qc.ld_step])
        writer.writerow(["ld_r2", qc.ld_r2])
        writer.writerow(["input_samples", len(samples)])
        writer.writerow(["input_sites", len(sites)])
        writer.writerow(["kept_samples", len(kept_samples)])
        writer.writerow(["kept_sites", len(kept_sites)])
        writer.writerow(["matrix_input_samples", matrix_stats.input_samples])
        writer.writerow(["matrix_input_sites", matrix_stats.input_sites])
        writer.writerow(["matrix_kept_samples", matrix_stats.kept_samples])
        writer.writerow(["matrix_kept_sites", matrix_stats.kept_sites])
        writer.writerow(["missingness_removed_sites", matrix_stats.missingness_removed_sites])
        writer.writerow(["monomorphic_removed_sites", matrix_stats.monomorphic_removed_sites])
        writer.writerow(["maf_removed_sites", matrix_stats.maf_removed_sites])
        writer.writerow(["sex_chr_removed_sites", matrix_stats.sex_chr_removed_sites])
        if ld_pruned_sites is not None:
            writer.writerow(["ld_pruned_sites", ld_pruned_sites])
    return out_path


def run_cohort_pca(config, samples: Iterable[str], *, force: bool = False) -> Dict[str, Path]:
    """Run resumable cohort PCA/MDS post-processing."""
    if config.data_type not in {"ancient", "modern", "auto"}:
        raise NotImplementedError("PCA stage supports data_type=ancient, modern, or auto")

    pca_sites = getattr(config.args, "pca_sites", None)
    if not pca_sites:
        raise ValueError("--pca-sites is required when --run-pca is specified")

    sample_list = sorted(set(samples))
    if len(sample_list) < 2:
        raise ValueError("PCA requires at least two successful samples")

    cohort_dir = config.results_dir / "cohort"
    pca_dir = cohort_dir / "pca"
    sites_table = cohort_dir / "pca_sites.tsv"
    resolved_data_type = infer_pca_data_type(config, sample_list, cohort_dir, force=force)
    if resolved_data_type == "ancient":
        raw_calls = cohort_dir / "pseudohaploid_raw_calls.tsv"
        matrix = cohort_dir / "pseudohaploid_matrix.tsv"
        filtered_matrix = cohort_dir / "pseudohaploid_matrix.filtered.tsv"
        ploidy = 1
        pseudo_haploid = True
    else:
        raw_calls = cohort_dir / "modern_genotype_calls.tsv"
        matrix = cohort_dir / "modern_genotype_matrix.tsv"
        filtered_matrix = cohort_dir / "modern_genotype_matrix.filtered.tsv"
        ploidy = 2
        pseudo_haploid = False
    qc = _pca_qc_from_args(config.args)
    engine = getattr(config.args, "pca_engine", "eigensoft")

    if _stage_done(sites_table, force):
        sites = read_pca_sites(sites_table)
    else:
        sites = read_pca_sites(Path(pca_sites))
        write_sites_table(sites, sites_table)

    if resolved_data_type == "ancient":
        sample_bams = {
            sample: _final_dedup_bam(config, sample)
            for sample in sample_list
        }
        if not _stage_done(raw_calls, force):
            extract_pseudohaploid_calls(
                sample_bams,
                sites,
                raw_calls,
                min_mapq=qc.min_mapq,
                min_baseq=qc.min_baseq,
                trim_ends=qc.trim_ends,
            )
        else:
            logger.info("再開: 既存pseudo-haploid raw callsを使用します: %s", raw_calls)

        if not _stage_done(matrix, force):
            build_pseudohaploid_matrix(raw_calls, sites, sample_list, matrix)
        else:
            logger.info("再開: 既存pseudo-haploid matrixを使用します: %s", matrix)
    else:
        sample_bams = {
            sample: _final_dedup_bam(config, sample)
            for sample in sample_list
        }
        if not _stage_done(raw_calls, force):
            extract_modern_diploid_calls(
                sample_bams,
                sites,
                raw_calls,
                min_mapq=qc.min_mapq,
                min_baseq=qc.min_baseq,
                trim_ends=qc.trim_ends,
            )
        else:
            logger.info("再開: 既存modern genotype callsを使用します: %s", raw_calls)

        if not _stage_done(matrix, force):
            build_modern_genotype_matrix(raw_calls, sites, sample_list, matrix)
        else:
            logger.info("再開: 既存modern genotype matrixを使用します: %s", matrix)

    if _stage_done(filtered_matrix, force):
        logger.info("再開: 既存filtered matrixを使用します: %s", filtered_matrix)
        matrix_stats = _filter_matrix_stats(
            matrix,
            qc.max_site_missing,
            qc.max_sample_missing,
            qc.min_maf,
            qc.exclude_sex_chr,
            sites,
            ploidy,
        )
    else:
        matrix_stats = _filter_matrix(
            matrix,
            filtered_matrix,
            qc.max_site_missing,
            qc.max_sample_missing,
            qc.min_maf,
            qc.exclude_sex_chr,
            sites,
            ploidy,
        )

    ld_pruned_sites: Optional[int] = None
    if engine == "eigensoft":
        eigenstrat_files, plink_files, pca_scores, pca_variance, mds, ld_pruned_sites = run_eigensoft_pca(
            filtered_matrix,
            sites,
            cohort_dir,
            qc,
            pseudo_haploid=pseudo_haploid,
            force=force,
        )
    else:
        eigenstrat_dir = cohort_dir / "eigenstrat"
        eigenstrat_files = (
            eigenstrat_dir / "cohort.geno",
            eigenstrat_dir / "cohort.snp",
            eigenstrat_dir / "cohort.ind",
        )
        if _stages_done(eigenstrat_files, force):
            logger.info("再開: 既存EIGENSTRAT互換出力を使用します: %s", eigenstrat_dir)
        else:
            eigenstrat_files = export_eigenstrat(filtered_matrix, sites, eigenstrat_dir)
        plink_dir = cohort_dir / "plink"
        plink_files = (plink_dir / "cohort.tped", plink_dir / "cohort.tfam")
        if _stages_done(plink_files, force):
            logger.info("再開: 既存PLINK text出力を使用します: %s", plink_dir)
        else:
            plink_files = export_plink_text(filtered_matrix, sites, plink_dir, pseudo_haploid=pseudo_haploid)
        pca_scores = pca_dir / "pca_scores.tsv"
        pca_variance = pca_dir / "pca_variance.tsv"
        mds = pca_dir / "mds.tsv"
        if _stages_done([pca_scores, pca_variance, mds], force):
            logger.info("再開: 既存Python PCA/MDS出力を使用します: %s", pca_dir)
        else:
            pca_scores, pca_variance, mds = run_pca_and_mds(filtered_matrix, pca_dir)
    qc_summary = _write_qc_summary(
        cohort_dir / "pca_qc_summary.tsv",
        sample_list,
        sites,
        filtered_matrix,
        qc,
        engine,
        matrix_stats,
        resolved_data_type,
        ploidy,
        config.data_type,
        ld_pruned_sites,
    )

    logger.info("cohort PCA/MDS stage が完了しました: %s", cohort_dir)
    return {
        "sites": sites_table,
        "raw_calls": raw_calls,
        "matrix": matrix,
        "filtered_matrix": filtered_matrix,
        "eigenstrat_geno": eigenstrat_files[0],
        "eigenstrat_snp": eigenstrat_files[1],
        "eigenstrat_ind": eigenstrat_files[2],
        "plink_tped": plink_files[0],
        "plink_tfam": plink_files[1],
        "pca_scores": pca_scores,
        "pca_variance": pca_variance,
        "mds": mds,
        "qc_summary": qc_summary,
    }
