import argparse
from pathlib import Path
import logging
from logging.handlers import RotatingFileHandler
import shutil
import subprocess
import sys
from typing import List, Optional, Union

from tqdm import tqdm

from tool_paths import CONVERTF_BIN, PICARD_JAR, PLINK_BIN, SMARTPCA_BIN

logger = logging.getLogger(__name__)

LOG_MAX_BYTES = 50 * 1024 * 1024
LOG_BACKUP_COUNT = 5

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
        # 2. fastq_dir/bam_dir が指定されている場合は project_accession を上書き
        self.fastq_dir: Optional[Path] = args.fastq_dir
        self.bam_dir: Optional[Path] = args.bam_dir
        if self.bam_dir:
            self.project_accession = self.bam_dir.name
        elif self.fastq_dir:
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
        # 6. 必要なディレクトリを自動生成
        for dir_path in [self.raw_data_dir, self.results_dir, self.logs_dir, self.temp_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

    def validate_environment(self) -> None:
        """
        パイプラインの実行に必要な外部ツール・ファイルの存在を一括検証する。
        不足があればまとめてエラーを表示し sys.exit(1) で終了する。
        """
        errors: List[str] = []

        required_tools = [
            "bwa", "samtools", "AdapterRemoval",
            "java", "gatk", "qualimap", "mapDamage",
        ]
        if getattr(self.args, "run_pca", False) and getattr(self.args, "pca_engine", "eigensoft") == "eigensoft":
            required_tools.extend([str(PLINK_BIN), str(CONVERTF_BIN), str(SMARTPCA_BIN)])
        for tool in required_tools:
            if not shutil.which(tool):
                errors.append(f"外部ツールが見つかりません: {tool}")

        if not PICARD_JAR.exists():
            errors.append(f"Picard jar が見つかりません: {PICARD_JAR}")

        bwa_index_suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
        missing_idx = [
            s for s in bwa_index_suffixes
            if not self.reference_genome.with_suffix(
                self.reference_genome.suffix + s
            ).exists()
        ]
        if missing_idx:
            errors.append(
                f"BWA インデックスが見つかりません (参照: {self.reference_genome}): "
                + ", ".join(missing_idx)
            )

        if errors:
            for msg in errors:
                logger.error(msg)
            logger.error(
                "環境検証に失敗しました。上記のツール/ファイルを確認してください。"
            )
            sys.exit(1)

        logger.info("環境検証に成功しました")


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
    # BAM ディレクトリ（ローカル、dedup BAM 入力）
    parser.add_argument(
        "--bam_dir",
        type=Path,
        default=None,
        help=(
            "既存の dedup BAM ファイルが存在するディレクトリへのパス。"
            "指定すると、FASTQ からの前処理をスキップして QC/VCF のみ実行します。"
        ),
    )
    parser.add_argument(
        "--bam_stage",
        type=str,
        default="dedup",
        choices=["dedup"],
        help="入力 BAM の処理段階。v1 では dedup のみ対応します。",
    )
    parser.add_argument(
        "--bam_pattern",
        type=str,
        default="*.bam",
        help="--bam_dir 配下で検出する BAM ファイルの glob パターン (デフォルト: *.bam)",
    )
    # データ種別
    parser.add_argument(
        "--data_type",
        type=str,
        default="ancient",
        help="データがmodernかancientか指定。auto は解析後にPCA分岐を推定します",
        choices=["modern", "ancient", "auto"],
    )
    # チェックポイントを無視して再実行
    parser.add_argument(
        "--force",
        action="store_true",
        help="完了済みサンプルのチェックポイント (.done) を無視して全サンプルを再実行する",
    )
    parser.add_argument(
        "--skip-vcf",
        action="store_true",
        help="sampleごとのGATK HaplotypeCaller VCF出力をスキップする。BAM QCとcohort PCA/MDSだけ実行したい場合に使う",
    )
    # リードグループ: ライブラリ名
    parser.add_argument(
        "--rg_library",
        type=str,
        default="unknown",
        help="Picard AddOrReplaceReadGroups の RGLB に設定するライブラリ名 (デフォルト: unknown)",
    )
    # リードグループ: シーケンシングセンター名
    parser.add_argument(
        "--rg_center",
        type=str,
        default="unknown",
        help="Picard AddOrReplaceReadGroups の RGCN に設定するシーケンシングセンター名 (デフォルト: unknown)",
    )
    # サンプル並列数
    parser.add_argument(
        "--parallel_samples",
        type=int,
        default=1,
        help=(
            "同時に解析するサンプル数。2 以上にするとサンプル単位で並列実行します。"
            "各サンプルが --threads 本のスレッドを使うため、合計スレッド数は "
            "parallel_samples × threads になります (デフォルト: 1)"
        ),
    )
    # HTTPS ダウンロード
    parser.add_argument(
        "--download-via-https",
        action="store_true",
        help="ena_download_https でダウンロードしてから解析を実行（再開対応・HTTPS）",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="端末上の進捗バー表示を無効化し、通常ログのみを出力する",
    )
    parser.add_argument(
        "--run-pca",
        action="store_true",
        help="全サンプル完了後に ancient DNA 向け cohort PCA/MDS stage を実行する",
    )
    parser.add_argument(
        "--pca-sites",
        type=Path,
        default=None,
        help="PCA/MDSで比較する共通SNPリスト (VCF または BED-like text)",
    )
    parser.add_argument(
        "--pca-engine",
        choices=["eigensoft", "python"],
        default="eigensoft",
        help="PCA実行エンジン。eigensoft は PLINK/CONVERTF/smartpca を使う (デフォルト: eigensoft)",
    )
    parser.add_argument(
        "--pca-min-mapq",
        type=int,
        default=30,
        help="pseudo-haploid抽出で使う最小mapping quality (デフォルト: 30)",
    )
    parser.add_argument(
        "--pca-min-baseq",
        type=int,
        default=30,
        help="pseudo-haploid抽出で使う最小base quality (デフォルト: 30)",
    )
    parser.add_argument(
        "--pca-trim-ends",
        type=int,
        default=2,
        help="pseudo-haploid抽出時にread両端から除外する塩基数 (デフォルト: 2)",
    )
    parser.add_argument(
        "--pca-max-sample-missing",
        type=float,
        default=0.9,
        help="PCA matrixで許容するsample欠損率の上限 (デフォルト: 0.9)",
    )
    parser.add_argument(
        "--pca-max-site-missing",
        type=float,
        default=0.9,
        help="PCA matrixで許容するSNP欠損率の上限 (デフォルト: 0.9)",
    )
    parser.add_argument(
        "--pca-min-maf",
        type=float,
        default=0.0,
        help="PCA matrixで残す最小minor allele frequency (デフォルト: 0.0)",
    )
    parser.add_argument(
        "--pca-exclude-sex-chr",
        action="store_true",
        help="PCA matrixから性染色体SNPを除外する",
    )
    parser.add_argument(
        "--pca-transversion-only",
        action="store_true",
        help="PCA matrixでtransversion SNPのみを残す。古DNA damage由来のC/T・G/A影響を抑える",
    )
    parser.add_argument(
        "--pca-ld-window",
        type=int,
        default=50,
        help="PLINK --indep-pairwise のwindow size (デフォルト: 50)",
    )
    parser.add_argument(
        "--pca-ld-step",
        type=int,
        default=5,
        help="PLINK --indep-pairwise のstep size (デフォルト: 5)",
    )
    parser.add_argument(
        "--pca-ld-r2",
        type=float,
        default=0.2,
        help="PLINK --indep-pairwise のr^2閾値 (デフォルト: 0.2)",
    )
    args = parser.parse_args()
    if args.bam_dir and args.fastq_dir:
        parser.error("--bam_dir と --fastq_dir は同時に指定できません")
    if args.bam_dir and args.download_via_https:
        parser.error("--bam_dir と --download-via-https は同時に指定できません")
    if args.run_pca and args.pca_sites is None:
        parser.error("--run-pca には --pca-sites が必要です")
    for name in ("pca_max_sample_missing", "pca_max_site_missing", "pca_min_maf", "pca_ld_r2"):
        value = getattr(args, name)
        if value < 0 or value > 1:
            parser.error("--%s は 0.0 から 1.0 の範囲で指定してください" % name.replace("_", "-"))
    for name in ("pca_min_mapq", "pca_min_baseq", "pca_trim_ends", "pca_ld_window", "pca_ld_step"):
        value = getattr(args, name)
        if value < 0:
            parser.error("--%s は 0 以上で指定してください" % name.replace("_", "-"))
    if args.pca_ld_window == 0 or args.pca_ld_step == 0:
        parser.error("--pca-ld-window と --pca-ld-step は 1 以上で指定してください")
    return args


class TqdmLoggingHandler(logging.StreamHandler):
    """tqdm の進捗バーを崩さずにコンソールへログを出す。"""

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
            tqdm.write(msg)
        except Exception:
            self.handleError(record)


def setup_logging(log_file: Path, use_tqdm: bool = False) -> logging.Logger:
    """
    パイプライン全体で使うロガーを初期化します。
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    stream_formatter = logging.Formatter(
        "%(asctime)s %(levelname)s: %(message)s"
    )
    file_formatter = logging.Formatter(
        "%(asctime)s %(levelname)s %(name)s: %(message)s"
    )

    has_stream_handler = any(
        isinstance(handler, logging.StreamHandler)
        and not isinstance(handler, logging.FileHandler)
        for handler in logger.handlers
    )
    if not has_stream_handler:
        stream_handler = TqdmLoggingHandler() if use_tqdm else logging.StreamHandler()
        stream_handler.setFormatter(stream_formatter)
        stream_handler.setLevel(logging.INFO)
        logger.addHandler(stream_handler)
    else:
        for handler in logger.handlers:
            if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                handler.setFormatter(stream_formatter)
                handler.setLevel(logging.INFO)

    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    file_exists = any(
        isinstance(handler, logging.FileHandler)
        and Path(handler.baseFilename) == log_file
        for handler in logger.handlers
    )
    if not file_exists:
        file_handler = RotatingFileHandler(
            log_file,
            maxBytes=LOG_MAX_BYTES,
            backupCount=LOG_BACKUP_COUNT,
        )
        file_handler.setFormatter(file_formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)
    else:
        for handler in logger.handlers:
            if isinstance(handler, logging.FileHandler) and Path(handler.baseFilename) == log_file:
                handler.setFormatter(file_formatter)
                handler.setLevel(logging.DEBUG)

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
