import argparse
from pathlib import Path
import logging
import re
import shlex
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
    return parser.parse_args()