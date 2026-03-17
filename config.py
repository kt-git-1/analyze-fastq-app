import argparse
from pathlib import Path
import logging
import re
import shlex
import shutil
import subprocess
from typing import Dict, List, Optional, Union


class PipelineConfig:
    """
    Pipeline configuration class は、ユーザーが指定した設定情報および、
    そこから導かれる設定値を一元管理するためのクラスです。
    """

    def __init__(self, args: argparse.Namespace) -> None:
        # 1. コマンドライン引数を取得
        self.args = args
        self.script_dir = Path(__file__).parent.resolve()
        self.base_dir: Path = args.base_dir
        # 2. fastq_dirが指定されている場合は project_accessionを上書き
        self.fastq_dir: Optional[Path] = args.fastq_dir
        if self.fastq_dir:
            self.project_accession = self.fastq_dir.name
        else:
            self.project_accession = args.project_accession
        # 3. ディレクトリ構造をセットアップ
        self.raw_data_dir: Path = self.base_dir / "raw_data" / self.project_accession
        self.results_dir: Path = self.base_dir / "results" / self.project_accession
        self.logs_dir: Path = self.base_dir / "logs"
        self.temp_dir: Path = self.base_dir / "temp"
        # 4. リファレンスゲノムのパスを決定
        self.reference_genome: Path = (
            args.reference_genome
            or (self.base_dir / "reference" / "equCab3.nochrUn.fa")
        )
        # 5. データタイプの保持
        self.data_type: str = getattr(args, "data_type", "ancient")
        # 7. 必要なディレクトリを自動生成
        for dir_path in [self.raw_data_dir, self.results_dir, self.logs_dir, self.temp_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)


def parse_args() -> argparse.Namespace:
    """
    コマンドライン引数を定義して解析するための関数。
    `main.py` で使用され、ユーザーが指定した設定を取得します。
    その結果はPipelineConfig()に渡されて、パイプライン全体の設定として使われます。

    Returns
    -------
    argparse.Namespace
        Parsed arguments namespace
    """
    script_dir = Path(__file__).parent.resolve()
    parser = argparse.ArgumentParser(
        description="NGSパイプラインの設定を管理するためのコマンドライン引数",
        fromfile_prefix_chars="@",
    )
    # プロジェクトアクセッション
    parser.add_argument(
        "--project_accession",
        default="PRJEB19970",
        help="ENAプロジェクトアクセッション (ローカルファイルを解析する場合は無視されます)",
    )
    # ベースディレクトリ
    parser.add_argument(
        "--base_dir",
        type=Path,
        default=script_dir / "data",
        help="データを保存するためのベースディレクトリ (raw, results, logs, temp)",
    )
    # 参照ゲノム
    parser.add_argument(
        "--reference_genome",
        type=Path,
        help="参照ゲノム FASTA ファイル",
    )
    # 並列ダウンロード数
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="並列ダウンロードのワーカー数（ローカル FASTQ モードの場合はほぼ使わない）",
    )
    # ダウンロードプロトコル
    parser.add_argument(
        "--download_protocol",
        type=str,
        choices=["ftp", "http"],
        default="http",
        help="FASTQ ファイルのダウンロードに使用するプロトコル。ftp または http を指定します (デフォルト: http)",
    )
    # ダウンロード失敗時の再試行回数
    parser.add_argument(
        "--max_retries",
        type=int,
        default=3,
        help="ダウンロードが失敗した場合の最大再試行回数 (デフォルト: 3)",
    )
    # 分析に使用するスレッド数
    parser.add_argument(
        "--threads",
        type=int,
        default=20,
        help="分析に使用するスレッド数（マッピング、バリアントコーリングなど）",
    )
    # Javaメモリ設定
    parser.add_argument(
        "--java_mem",
        default="10g",
        help="GATKやその他のJavaツールに使用するJavaメモリ設定",
    )
    # FASTQ ディレクトリ（ローカル）
    parser.add_argument(
        "--fastq_dir",
        type=Path,
        default=None,
        help=(
            "事前にダウンロード済みの FASTQ ファイルが存在するディレクトリへのパス。 "
            "指定すると、パイプラインは ENA からのダウンロードをスキップし、 "
            "このディレクトリからペアエンドの FASTQ ファイルを直接読み込みます。"
        ),
    )
    # データ種別
    parser.add_argument(
        "--data_type",
        type=str,
        default="ancient",
        help="データがmodernかancientか指定",
        choices=["modern", "ancient"],
    )
    # HTTPS ダウンロード
    parser.add_argument(
        "--download-via-https",
        action="store_true",
        help="ena_download_https でダウンロードしてから解析を実行（再開対応・HTTPS）",
    )
    return parser.parse_args()


def setup_logging(log_file: Path) -> logging.Logger:
    """
    パイプライン全体で使うロガーを初期化します。
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(name)s: %(message)s"
    )

    has_stream_handler = any(
        isinstance(handler, logging.StreamHandler)
        and not isinstance(handler, logging.FileHandler)
        for handler in logger.handlers
    )
    if not has_stream_handler:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    file_exists = any(
        isinstance(handler, logging.FileHandler)
        and Path(handler.baseFilename) == log_file
        for handler in logger.handlers
    )
    if not file_exists:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


def cleanup_intermediate_file(file_path: Union[str, Path], logger: logging.Logger) -> None:
    """
    中間ファイルを削除します。存在しない場合は何もしません。
    """
    path = Path(file_path)
    if not path.exists():
        return

    try:
        path.unlink()
        logger.info("中間ファイルを削除しました: %s", path)
    except OSError as exc:
        logger.warning("中間ファイルを削除できませんでした: %s (%s)", path, exc)


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

        # sample ごとにサブディレクトリが切られている場合は最上位ディレクトリ名を優先。
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