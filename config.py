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

    def __init__(self, args: argparse.Namespace):
        # 1. コマンドライン引数を取得
        self.args = args
        self.script_dir = Path(__file__).parent.resolve()
        self.base_dir: Path = args.base_dir

        # 2. fastq_dirが指定されている場合は project_accessionを上書き
        self.fastq_dir: Optional[Path] = args.fastq_dir
        if self.fastq_dir:
            self.project_accession = self.fastq_dir.name
        else:
            self.project_accession: str = args.project_accession

        # 3. ディレクトリ構造をセットアップ
        # ローカルファイルを解析する場合、raw_data_dirは使用されないが、下流のモジュールとの互換性を維持するために保持されています。
        # 結果とログは、データソースに関係なく、base_dirに対して相対的に書き込まれます。
        self.raw_data_dir: Path = self.base_dir / "raw_data" / self.project_accession
        self.results_dir: Path = self.base_dir / "results" / self.project_accession
        self.logs_dir: Path = self.base_dir / "logs"
        self.temp_dir: Path = self.base_dir / "temp"

        # 4. リファレンスゲノムのパスを決定
        # デフォルトでは `equCab3.nochrUn.fa` を設定
        self.reference_genome: Path = (
            args.reference_genome
            or (self.base_dir / "reference" / "equCab3.nochrUn.fa")
        )

        # 5. ローカルの FASTQ ディレクトリ（任意）を保存
        # ローカルファイルを解析する場合、ENAからのダウンロードはスキップされ、代わりにこのディレクトリを探索してR1/R2ペアを検出します。
        self.fastq_dir: Optional[Path] = args.fastq_dir

        # 6. 必要なディレクトリを自動生成
        # ローカルファイルを解析する場合、raw_data_dirは使用されないが、下流のモジュールとの互換性を維持するために保持されています。
        for dir_path in [
            self.raw_data_dir,
            self.results_dir,
            self.logs_dir,
            self.temp_dir,
        ]:
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
    # 1. スクリプトのディレクトリを取得
    script_dir = Path(__file__).parent.resolve()
    # 2. 引数パーサーを作成
    parser = argparse.ArgumentParser(
        description="NGSパイプラインの設定を管理するためのコマンドライン引数",
        fromfile_prefix_chars="@",
    )
    # 3. プロジェクトアクセションを定義
    parser.add_argument(
        "--project_accession",
        default="PRJEB19970",
        help="ENAプロジェクトアクセション (ローカルファイルを解析する場合は無視されます)",
    )
    # 4. ベースディレクトリを定義
    parser.add_argument(
        "--base_dir",
        type=Path,
        default=script_dir / "data",
        help="データを保存するためのベースディレクトリ (raw, results, logs, temp)",
    )
    # 5. リファレンスゲノムのパスを定義
    parser.add_argument(
        "--reference_genome",
        type=Path,
        help="参照ゲノム FASTA ファイル",
    )
    # 6. ダウンロードの並列数を定義
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="並列ダウンロードのワーカー数（ローカル FASTQ モードの場合はほぼ使わない）",
    )
    # 7. スレッド数を定義
    parser.add_argument(
        "--threads",
        type=int,
        default=20,
        help="分析に使用するスレッド数（マッピング、バリアントコーリングなど）",
    )
    # 8. Javaメモリ設定を定義
    parser.add_argument(
        "--java_mem",
        default="10g",
        help="GATKやその他のJavaツールに使用するJavaメモリ設定",
    )
    # 9. ローカルの FASTQ ディレクトリを定義
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
    return parser.parse_args()


def setup_logging(log_file: Optional[Path] = None) -> logging.Logger:
    """
    パイプラインのロギングを初期化するための関数。
    ログは標準出力と、指定された場合はログファイルに書き込まれます。
    ログファイルは、パイプラインの実行中に生成されます。

    Parameters
    ----------
    log_file : Optional[Path], optional
        ログファイルへのパス（指定されない場合は標準出力のみ）

    Returns
    -------
    logging.Logger
        設定されたLoggerインスタンス
    """
    # 1. ハンドラーを初期化
    handlers = []
    # 2. コンソールハンドラーを追加
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    handlers.append(console_handler)
    # 3. ログファイルハンドラーを追加（指定された場合）
    if log_file:
        file_handler = logging.FileHandler(str(log_file))
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        )
        handlers.append(file_handler)
    # 4. ロギングを初期化
    logging.basicConfig(level=logging.INFO, handlers=handlers)
    # 5. Loggerインスタンスを取得
    logger = logging.getLogger(__name__)
    if log_file:
        logger.info(f"ログをファイルに出力しています: {log_file}")
    return logger


def cleanup_intermediate_file(path: Path, logger: logging.Logger) -> None:
    """
    中間ファイルを削除するための関数。
    削除に失敗した場合は、パイプラインを中断せずに警告ログを出力します。

    Parameters
    ----------
    path : Path
        削除するファイルのパス
    logger : logging.Logger
        ステータスレポートのためのLogger
    """
    # 1. ファイルが存在しない場合は何もしない
    if not path:
        return
    try:
        # 2. ファイルが存在する場合は削除
        if path.exists():
            path.unlink()
            logger.info(f"中間ファイルを削除しました: {path}")
    except Exception as exc:
        logger.warning(f"中間ファイルを削除できませんでした: {path}: {exc}")


def parse_fastq_general(directory: Path) -> Dict[str, Dict[str, List[Path]]]:
    """
    ZYJ2系, Filgen系, シングルトン系に対応して、
    サンプルID × R1/R2 ごとに FASTQ/FQ をまとめる。

    戻り値:
      {
        sample_id: {
          "R1": [Path(...), ...],
          "R2": [Path(...), ...],  # 無ければキー無し or 空リスト
        },
        ...
      }
    """

    # 例1: ZYJ2_S1_L005_R1_001.fastq.gz
    illumina_pat = re.compile(
        r"^(.+)_L\d{3}_R([12])_001\.fastq\.gz$"
    )

    # 例2: J01_001_NDSW42296_H7VHWDSXX_L1_1.fq.gz
    filgen_pat = re.compile(
        r"^(.+)_L\d+_([12])\.f(?:ast)?q\.gz$"
    )

    # 例3: シングルトン（レーン情報やR1/R2がないもの）
    # BER01_A__BER01_A_SCY2.2_user_TGACGT__... .fastq.gz など
    singleton_pat = re.compile(
        r"^(.+)\.f(?:ast)?q\.gz$"
    )

    sample_to_reads: Dict[str, Dict[str, List[Path]]] = {}

    # サブディレクトリも含めて探索
    for fastq in directory.rglob("*.f*q.gz"):
        name = fastq.name

        sample_id = None
        read_num = None  # "1" / "2"

        # 1. Illumina パターン
        m = illumina_pat.match(name)
        if m:
            sample_id, read_num = m.group(1), m.group(2)
        else:
            # 2. Filgen パターン
            m = filgen_pat.match(name)
            if m:
                # sample_id をファイル名ベースにする場合
                sample_id, read_num = m.group(1), m.group(2)
                # フォルダ名をサンプルIDにしたいなら ↓ を使う
                # sample_id = fastq.parent.name
            else:
                # 3. シングルトン → R1 とみなす
                m = singleton_pat.match(name)
                if m:
                    sample_id = m.group(1)  # 拡張子を除いた全部をサンプルIDに
                    read_num = "1"
                else:
                    # どのパターンにも合わなければスキップ
                    continue

        reads = sample_to_reads.setdefault(sample_id, {})
        key = f"R{read_num}"           # "R1" or "R2"
        reads.setdefault(key, []).append(fastq)

    # 後段処理のため、ファイル名でソートしておく
    for sample_id, reads in sample_to_reads.items():
        for key in ("R1", "R2"):
            if key in reads:
                reads[key].sort(key=lambda p: p.name)

    return sample_to_reads



def merge_lanes_by_cat(
    sample_to_reads: Dict[str, Dict[str, List[Path]]],
    out_dir: Path,
    logger: Optional[object] = None,
) -> Dict[str, List[Path]]:
    """
    サンプル × R1/R2 ごとに複数レーンがあれば cat でマージする。
    1ファイルしかない場合は元ファイルをそのまま返す（ファイル名に merged を付けない）。

    Returns
    -------
    Dict[str, List[Path]]
        {
          sample_id: [R1_path, R2_path]  # R2 がなければ [R1_path]
        }
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    merged: Dict[str, List[Path]] = {}

    for sample_id, reads in sample_to_reads.items():
        result_files: List[Path] = []

        for read_key in ("R1", "R2"):
            files = reads.get(read_key)
            if not files:
                continue

            files = sorted(files, key=lambda p: p.name)

            # ★ 1ファイル → cat しない・merged名にしない
            if len(files) == 1:
                if logger:
                    logger.info(f"{sample_id} {read_key}: single FASTQ → no cat")
                result_files.append(files[0])
                continue

            # ★ 複数ファイル → cat する
            suffix = "".join(files[0].suffixes)  # .fastq.gz / .fq.gz
            merged_path = out_dir / f"{sample_id}_{read_key}_merged{suffix}"

            file_list = " ".join(shlex.quote(str(p)) for p in files)
            cmd = f"cat {file_list} > {shlex.quote(str(merged_path))}"

            if logger:
                logger.info(f"{sample_id} {read_key}: merging {len(files)} files → {merged_path}")
                logger.info(f"Command: {cmd}")

            subprocess.run(["bash", "-lc", cmd], check=True)

            result_files.append(merged_path)

        if result_files:
            merged[sample_id] = result_files

    return merged