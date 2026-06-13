import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PCASite:
    chrom: str
    pos: int
    site_id: str
    ref: str = ""
    alt: str = ""


def _is_nonempty_file(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _stage_done(path: Path, force: bool) -> bool:
    return (not force) and _is_nonempty_file(path)


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

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "site_id", "chrom", "pos", "ref", "alt", "allele", "baseq", "mapq"])
        for sample, bam_path in sorted(sample_bams.items()):
            if not _is_nonempty_file(bam_path):
                logger.warning("PCA allele extraction skipped missing BAM: %s (%s)", sample, bam_path)
                continue
            logger.info("pseudo-haploid allele を抽出します: %s", sample)
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for site in sites:
                    allele = ""
                    baseq = 0
                    mapq = 0
                    for column in bam.pileup(
                        site.chrom,
                        site.pos - 1,
                        site.pos,
                        truncate=True,
                        stepper="samtools",
                        min_mapping_quality=min_mapq,
                    ):
                        if column.reference_pos != site.pos - 1:
                            continue
                        allele, baseq, mapq = _choose_pileup_allele(
                            column,
                            site,
                            min_baseq,
                            trim_ends,
                        )
                        break
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


def filter_matrix(
    matrix_path: Path,
    out_path: Path,
    *,
    max_site_missing: float = 0.9,
    max_sample_missing: float = 0.9,
) -> Path:
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
    variable = []
    for col_idx in range(data.shape[1]):
        values = data[:, col_idx]
        observed = values[~np.isnan(values)]
        variable.append(len(set(observed.tolist())) > 1)
    keep_sites = (site_missing <= max_site_missing) & np.array(variable, dtype=bool)
    data = data[:, keep_sites]
    kept_site_ids = [site_id for site_id, keep in zip(header[1:], keep_sites) if keep]

    if data.size == 0 or not kept_site_ids:
        raise ValueError("No variable SNP sites remain after PCA filtering")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + kept_site_ids)
        for row, values in zip(kept_rows, data):
            writer.writerow([row[0]] + ["" if np.isnan(value) else str(int(value)) for value in values])
    return out_path


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


def export_plink_text(matrix_path: Path, sites: List[PCASite], out_dir: Path) -> Tuple[Path, Path]:
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
                else:
                    row.extend([site.alt or "B", site.alt or "B"])
            writer.writerow(row)
    with tfam_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=" ")
        for sample in samples:
            writer.writerow([sample, sample, "0", "0", "0", "-9"])
    return tped_path, tfam_path


def _write_qc_summary(
    out_path: Path,
    samples: List[str],
    sites: List[PCASite],
    filtered_matrix: Path,
) -> Path:
    kept_samples, kept_sites, _data = _load_numeric_matrix(filtered_matrix)
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["input_samples", len(samples)])
        writer.writerow(["input_sites", len(sites)])
        writer.writerow(["kept_samples", len(kept_samples)])
        writer.writerow(["kept_sites", len(kept_sites)])
    return out_path


def run_cohort_pca(config, samples: Iterable[str], *, force: bool = False) -> Dict[str, Path]:
    """Run resumable ancient-DNA-oriented cohort PCA/MDS post-processing."""
    if config.data_type != "ancient":
        raise NotImplementedError("PCA stage currently supports data_type=ancient only")

    pca_sites = getattr(config.args, "pca_sites", None)
    if not pca_sites:
        raise ValueError("--pca-sites is required when --run-pca is specified")

    sample_list = sorted(set(samples))
    if len(sample_list) < 2:
        raise ValueError("PCA requires at least two successful samples")

    cohort_dir = config.results_dir / "cohort"
    pca_dir = cohort_dir / "pca"
    sites_table = cohort_dir / "pca_sites.tsv"
    raw_calls = cohort_dir / "pseudohaploid_raw_calls.tsv"
    matrix = cohort_dir / "pseudohaploid_matrix.tsv"
    filtered_matrix = cohort_dir / "pseudohaploid_matrix.filtered.tsv"

    if _stage_done(sites_table, force):
        sites = read_pca_sites(sites_table)
    else:
        sites = read_pca_sites(Path(pca_sites))
        write_sites_table(sites, sites_table)

    sample_bams = {
        sample: config.results_dir / sample / "dedup" / "%s.dedup.sorted.bam" % sample
        for sample in sample_list
    }

    if not _stage_done(raw_calls, force):
        extract_pseudohaploid_calls(sample_bams, sites, raw_calls)
    else:
        logger.info("再開: 既存pseudo-haploid raw callsを使用します: %s", raw_calls)

    if not _stage_done(matrix, force):
        build_pseudohaploid_matrix(raw_calls, sites, sample_list, matrix)
    else:
        logger.info("再開: 既存pseudo-haploid matrixを使用します: %s", matrix)

    if not _stage_done(filtered_matrix, force):
        filter_matrix(matrix, filtered_matrix)
    else:
        logger.info("再開: 既存filtered matrixを使用します: %s", filtered_matrix)

    eigenstrat_files = export_eigenstrat(filtered_matrix, sites, cohort_dir / "eigenstrat")
    plink_files = export_plink_text(filtered_matrix, sites, cohort_dir / "plink")
    pca_scores, pca_variance, mds = run_pca_and_mds(filtered_matrix, pca_dir)
    qc_summary = _write_qc_summary(cohort_dir / "pca_qc_summary.tsv", sample_list, sites, filtered_matrix)

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
