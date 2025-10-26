import re
import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, List

import requests
from tqdm import tqdm

from config import PipelineConfig, parse_args, setup_logging
from modules.ena_downloader import ENADownloader
from modules.bwa_mapper import BWAMapper
from modules.softclipper import SoftClipper
from modules.bam_processor import BAMProcessor
from modules.analyzers import (
    MapDamageAnalyzer,
    QualimapAnalyzer,
    HaplotypeCaller,
)


def _cleanup_intermediate_file(path: Path, logger: logging.Logger) -> None:
    """
    Remove an intermediate file if it exists.  Failures to remove
    intermediate files should not halt the pipeline, hence any
    exceptions are caught and logged as warnings.

    Parameters
    ----------
    path : Path
        Path to the file to remove.
    logger : logging.Logger
        Logger for status reporting.
    """
    if not path:
        return
    try:
        if path.exists():
            path.unlink()
            logger.info(f"Removed intermediate file: {path}")
    except Exception as exc:
        logger.warning(f"Failed to remove intermediate file {path}: {exc}")


def parse_local_fastq(directory: Path) -> Dict[str, List[Path]]:
    """
    Discover paired FASTQ files within a directory based on Illumina
    naming conventions.  Files are grouped by the common prefix
    preceding the ``_R1_001.fastq.gz`` or ``_R2_001.fastq.gz`` suffix.

    The grouping logic assumes that paired reads share the same lane
    (``L###``) and sample identifiers, differing only by the ``R1`` or
    ``R2`` read indicator.  For example, ``PE-AncientHorses-01_S1_L001_R1_001.fastq.gz``
    and ``PE-AncientHorses-01_S1_L001_R2_001.fastq.gz`` constitute a
    pair, while ``PE-AncientHorses-01_S1_L002_R1_001.fastq.gz`` and
    ``PE-AncientHorses-01_S1_L002_R2_001.fastq.gz`` form another pair.

    Parameters
    ----------
    directory : Path
        Directory containing FASTQ files.

    Returns
    -------
    Dict[str, List[Path]]
        Mapping from a sample identifier (derived from the filename
        prefix) to a list of one or two FASTQ files. The order within
        each list is ``[R1, R2]`` when both reads are present, or
        ``[R1]`` when only a single‑end read is available.
    """
    sample_to_files: Dict[str, Dict[str, Path]] = {}
    pattern = re.compile(r"(.+)_R([12])_001\.fastq\.gz$")
    for fastq in directory.glob("*.fastq.gz"):
        match = pattern.match(fastq.name)
        if not match:
            # Skip files that do not follow the expected naming scheme
            continue
        prefix, read_number = match.groups()
        entry = sample_to_files.setdefault(prefix, {})
        if read_number == "1":
            entry["R1"] = fastq
        elif read_number == "2":
            entry["R2"] = fastq
    result: Dict[str, List[Path]] = {}
    for prefix, reads in sample_to_files.items():
        # Ensure deterministic ordering: R1 first, then R2 if present
        ordered: List[Path] = []
        if "R1" in reads:
            ordered.append(reads["R1"])
        if "R2" in reads:
            ordered.append(reads["R2"])
        result[prefix] = ordered
    return result


def main() -> None:
    """
    Entry point for the complete NGS analysis pipeline.  Depending on
    user input, this function either downloads reads from ENA or
    processes pre‑downloaded FASTQ files from a local directory.  For
    each sample or lane, the pipeline performs adapter trimming,
    mapping, soft clipping, BAM processing, quality control analyses,
    and variant calling.  Intermediate files are cleaned up after
    downstream analyses complete successfully.
    """
    # Initialise configuration and logging
    args = parse_args()
    config = PipelineConfig(args)
    log_file = config.logs_dir / f"pipeline_{config.project_accession}.log"
    logger = setup_logging(log_file=log_file)

    logger.info(
        f"Starting complete pipeline for project {config.project_accession}"
    )

    # Verify presence of the reference genome
    if not config.reference_genome.exists():
        logger.error(f"Reference genome not found: {config.reference_genome}")
        sys.exit(1)

    # Build reference index if missing
    fai = config.reference_genome.with_suffix(".fai")
    if not fai.exists():
        logger.info("Creating reference genome index")
        subprocess.run(["samtools", "faidx", str(config.reference_genome)], check=True)

    # Instantiate modules
    ena_downloader = ENADownloader(config)
    bwa_mapper = BWAMapper(config)
    softclipper = SoftClipper(config)
    bam_processor = BAMProcessor(config)
    mapdamage_analyzer = MapDamageAnalyzer(config)
    qualimap_analyzer = QualimapAnalyzer(config)
    haplotypecaller = HaplotypeCaller(config)

    session = requests.Session()

    # Determine data source: either ENA or local FASTQ directory
    if config.fastq_dir:
        # Local FASTQ processing
        if not config.fastq_dir.exists() or not config.fastq_dir.is_dir():
            logger.error(
                f"Specified FASTQ directory does not exist or is not a directory: {config.fastq_dir}"
            )
            sys.exit(1)
        logger.info(f"Using local FASTQ files from: {config.fastq_dir}")
        sample_to_fastqs = parse_local_fastq(config.fastq_dir)
    else:
        # Download from ENA as before
        logger.info("Step 1: Downloading data from ENA")
        response_data = ena_downloader.get_api_response(
            config.project_accession, session
        )
        sample_to_ftp_urls = ena_downloader.parse_response_data(response_data)
        sample_to_fastqs = {}
        # Download FASTQ files for each sample
        for sample_acc, ftp_urls in sample_to_ftp_urls.items():
            logger.info(f"Processing sample: {sample_acc}")
            fastq_files = ena_downloader.download_sample_data(sample_acc, ftp_urls)
            if not fastq_files:
                logger.warning(f"No FASTQ files downloaded for {sample_acc}")
                continue
            # Sort to ensure deterministic ordering; R1 then R2
            fastq_files_sorted: List[Path] = sorted(
                fastq_files,
                key=lambda p: ("R2" in p.name, p.name),
            )
            sample_to_fastqs[sample_acc] = fastq_files_sorted

    # Iterate over samples (either downloaded or local)
    for sample_acc, fastq_files in tqdm(
        sample_to_fastqs.items(), desc="progress", unit="sample"
    ):
        logger.info(f"Running analysis for sample/lane: {sample_acc}")
        # Step 2: BWA mapping (paired or single)
        bam_file = bwa_mapper.run_mapping_pipeline(sample_acc, fastq_files)
        if not bam_file:
            logger.error(f"Mapping failed for {sample_acc}")
            continue
        # Step 3: Soft clipping
        softclipped_bam = softclipper.run_softclipping(sample_acc, bam_file)
        if not softclipped_bam:
            logger.error(f"Soft clipping failed for {sample_acc}")
            continue
        _cleanup_intermediate_file(bam_file, logger)
        # Step 4: BAM processing (sort, dedup, index)
        dedup_bam = bam_processor.run_bam_processing(sample_acc, softclipped_bam)
        if not dedup_bam:
            logger.error(f"BAM processing failed for {sample_acc}")
            continue
        # Step 5: mapDamage analysis
        mapdamage_result = mapdamage_analyzer.run_mapdamage(sample_acc, softclipped_bam)
        if not mapdamage_result:
            logger.error(f"MapDamage analysis failed for {sample_acc}")
            continue
        # Clean up BWA intermediate files created in bam_files/ directory
        bam_dir = dedup_bam.parent.parent / "bam_files"
        if bam_dir.exists():
            for pattern in ("*.bam", "*.bai", "*.truncated"):
                for intermediate in bam_dir.glob(pattern):
                    _cleanup_intermediate_file(intermediate, logger)
        # Step 6: Qualimap analysis
        qualimap_result = qualimap_analyzer.run_qualimap(sample_acc, dedup_bam)
        if not qualimap_result:
            logger.error(f"Qualimap analysis failed for {sample_acc}")
            continue
        # Step 7: HaplotypeCaller
        vcf_file = haplotypecaller.run_haplotypecaller(sample_acc, dedup_bam)
        if not vcf_file:
            logger.error(f"HaplotypeCaller failed for {sample_acc}")
            continue
        # Cleanup deduplicated BAM and its index after analyses are complete
        dedup_bam_index = Path(str(dedup_bam) + ".bai")
        _cleanup_intermediate_file(dedup_bam_index, logger)
        _cleanup_intermediate_file(dedup_bam, logger)
        _cleanup_intermediate_file(softclipped_bam, logger)
        logger.info(f"Completed processing for {sample_acc}")
    logger.info("Pipeline completed successfully!")


if __name__ == "__main__":
    main()