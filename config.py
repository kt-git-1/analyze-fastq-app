import argparse
from pathlib import Path
import logging
import shutil
import subprocess
import sys
from typing import List, Optional, Union

logger = logging.getLogger(__name__)


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
        for tool in required_tools:
            if not shutil.which(tool):
                errors.append(f"外部ツールが見つかりません: {tool}")

        picard_jar: Path = Path(self.args.picard_jar)
        if not picard_jar.exists():
            errors.append(f"Picard jar が見つかりません: {picard_jar}")

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
    # データ種別
    parser.add_argument(
        "--data_type",
        type=str,
        default="ancient",
        help="データがmodernかancientか指定",
        choices=["modern", "ancient"],
    )
    # チェックポイントを無視して再実行
    parser.add_argument(
        "--force",
        action="store_true",
        help="完了済みサンプルのチェックポイント (.done) を無視して全サンプルを再実行する",
    )
    # Picard jar パス
    parser.add_argument(
        "--picard_jar",
        type=Path,
        default=Path("/usr/local/bin/picard.jar"),
        help="Picard jar ファイルのパス (デフォルト: /usr/local/bin/picard.jar)",
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