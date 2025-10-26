import argparse
from pathlib import Path
import logging
from typing import Optional


class PipelineConfig:
    """
    Pipeline configuration class used to centralize all user-provided and
    derived settings. Previously this project only supported pulling read
    data directly from the European Nucleotide Archive (ENA). To support
    local analysis of pre‑downloaded FASTQ files, a new argument
    ``--fastq_dir`` has been introduced. When specified, the pipeline
    will bypass any remote download steps and instead parse paired FASTQ
    files from the provided directory. See :func:`parse_local_fastq` in
    ``main.py`` for details on the file grouping logic.
    """

    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.script_dir = Path(__file__).parent.resolve()
        self.base_dir: Path = args.base_dir
        self.project_accession: str = args.project_accession

        # Directory structure
        # When analysing local files, raw_data_dir may not be used, but it
        # remains available to maintain compatibility with downstream
        # modules (e.g. for logging or archiving).  Results and logs are
        # still written relative to base_dir irrespective of the data source.
        self.raw_data_dir: Path = self.base_dir / "raw_data" / self.project_accession
        self.results_dir: Path = self.base_dir / "results" / self.project_accession
        self.logs_dir: Path = self.base_dir / "logs"
        self.temp_dir: Path = self.base_dir / "temp"

        # Reference genome path. If not provided on the command line this
        # defaults to ``equCab3.fa`` as before.
        self.reference_genome: Path = (
            args.reference_genome
            or (self.base_dir / "reference" / "equCab3.fa")
        )

        # Optional directory containing pre‑downloaded FASTQ files. When
        # present, ``main.py`` will ignore ``ENA`` downloads and instead
        # iterate through this directory to discover R1/R2 pairs.
        self.fastq_dir: Optional[Path] = args.fastq_dir

        # Create required directories.  ``raw_data_dir`` may already exist
        # if local data lives elsewhere.  ``exist_ok=True`` prevents errors.
        for dir_path in [
            self.raw_data_dir,
            self.results_dir,
            self.logs_dir,
            self.temp_dir,
        ]:
            dir_path.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments for the pipeline.  This function was
    extended with an additional ``--fastq_dir`` argument to support
    processing of pre‑downloaded FASTQ files. When this option is
    supplied, the pipeline will skip the ENA download step and
    instead look for files in the provided directory.

    Returns
    -------
    argparse.Namespace
        Parsed arguments namespace
    """
    script_dir = Path(__file__).parent.resolve()
    parser = argparse.ArgumentParser(
        description="Complete ENA download and analysis pipeline",
        fromfile_prefix_chars="@",
    )
    parser.add_argument(
        "--project_accession",
        default="PRJEB19970",
        help="ENA project accession (ignored when --fastq_dir is provided)",
    )
    parser.add_argument(
        "--base_dir",
        type=Path,
        default=script_dir / "data",
        help="Base directory for all data (raw, results, logs, temp)",
    )
    parser.add_argument(
        "--reference_genome",
        type=Path,
        help="Reference genome FASTA file",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel download workers",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=20,
        help="Number of threads for analysis (mapping, variant calling, etc.)",
    )
    parser.add_argument(
        "--java_mem",
        default="10g",
        help="Java memory setting for GATK and other Java tools",
    )
    parser.add_argument(
        "--fastq_dir",
        type=Path,
        default=None,
        help=(
            "Path to a directory containing pre‑downloaded FASTQ files. "
            "If provided, the pipeline will skip ENA downloads and read "
            "paired FASTQ files directly from this directory."
        ),
    )
    return parser.parse_args()


def setup_logging(log_file: Optional[Path] = None) -> logging.Logger:
    """
    Initialize logging for the pipeline.  Logs are written both to
    standard output and, if provided, to a log file.  The presence of
    this function remains unchanged from the original version.

    Parameters
    ----------
    log_file : Optional[Path], optional
        Path to a log file, by default None

    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    handlers = []
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    handlers.append(console_handler)
    # File handler (if a log file is specified)
    if log_file:
        file_handler = logging.FileHandler(str(log_file))
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        )
        handlers.append(file_handler)
    logging.basicConfig(level=logging.INFO, handlers=handlers)
    logger = logging.getLogger(__name__)
    if log_file:
        logger.info(f"ログをファイルに出力しています: {log_file}")
    return logger