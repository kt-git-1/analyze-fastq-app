import hashlib
import requests
import logging
import os
import sys
import time
from ftplib import FTP
from urllib.parse import urlparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# Set up a module-level logger. 既存のロガーを使う場合は上書きされないようにします。
logger = logging.getLogger(__name__)


def verify_file_md5(path: Path, expected: str) -> bool:
    """ファイルの MD5 ハッシュを計算し、期待値と比較する。"""
    md5 = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            md5.update(chunk)
    actual = md5.hexdigest()
    if actual != expected:
        logger.warning(
            "MD5 mismatch: %s (expected=%s, actual=%s)", path.name, expected, actual
        )
        return False
    return True


class ENADownloader:
    """
    European Nucleotide Archive (ENA) から FASTQ ファイルをダウンロードするためのクラス。

    デフォルトでは FTP URL を HTTP に変換してダウンロードを試み、失敗した場合はリトライを行います。
    コマンドライン引数で `--download_protocol` と `--max_retries` を設定すると挙動を変更できます。
    """

    def __init__(self, config) -> None:
        """
        Initialize downloader with pipeline configuration.

        Parameters
        ----------
        config : PipelineConfig
            パイプライン全体の設定オブジェクト。
        """
        self.config = config
        # Use command-line arguments if they exist, otherwise default values.
        self.protocol: str = getattr(config.args, "download_protocol", "http")
        self.max_retries: int = getattr(config.args, "max_retries", 3)

    def get_api_response(self, project_accession: str, session: requests.Session) -> str:
        """
        ENA Portal API からプロジェクト単位でランごとのメタデータを取得します。

        Parameters
        ----------
        project_accession : str
            対象のプロジェクトアクセッション（例: PRJEB19970）。
        session : requests.Session
            再利用可能な HTTP セッション。

        Returns
        -------
        str
            タブ区切りのファイルレポート。
        """
        url = (
            "https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={project_accession}"
            "&result=read_run"
            "&fields=sample_accession,submitted_ftp,fastq_md5"
            "&format=tsv"
        )
        try:
            response = session.get(url, timeout=10)
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as exc:
            # API 呼び出しに失敗した場合はエラーログを出力して終了
            logger.error(f"プロジェクト {project_accession} のデータ取得に失敗しました: {exc}")
            sys.exit(1)

    def parse_response_data(self, response_data: str) -> dict:
        """
        API レスポンスを解析してサンプルごとの FTP URL のリストを作成します。

        Parameters
        ----------
        response_data : str
            API から受け取ったタブ区切りテキスト。

        Returns
        -------
        dict
            {sample_accession: [ftp_url, ...], ...} の形式の辞書。
        """
        sample_to_ftp_urls: dict = {}
        # 先頭行はヘッダーなので除外
        lines = response_data.strip().split("\n")[1:]
        for line in lines:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            sample_acc, ftp_urls = parts[0], parts[1]
            for ftp_url in ftp_urls.split(";"):
                ftp_url = ftp_url.strip()
                if not ftp_url:
                    continue
                sample_to_ftp_urls.setdefault(sample_acc, []).append(ftp_url)
        return sample_to_ftp_urls

    def parse_response_with_checksums(self, response_data: str) -> dict:
        """
        API レスポンスを解析してサンプルごとの (URL, MD5) タプルのリストを作成します。

        Returns
        -------
        dict
            {sample_accession: [(ftp_url, md5_or_None), ...], ...}
        """
        sample_to_files: dict = {}
        lines = response_data.strip().split("\n")
        if len(lines) < 2:
            return sample_to_files

        header = lines[0].split("\t")
        try:
            sample_idx = header.index("sample_accession")
            url_idx = header.index("submitted_ftp")
        except ValueError:
            logger.error("API レスポンスのヘッダーが不正です: %s", header)
            return sample_to_files

        md5_idx = header.index("fastq_md5") if "fastq_md5" in header else None

        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) <= max(sample_idx, url_idx):
                continue

            sample_acc = parts[sample_idx]
            ftp_urls = parts[url_idx].split(";")
            md5s = (
                parts[md5_idx].split(";")
                if md5_idx is not None and len(parts) > md5_idx
                else []
            )

            for i, ftp_url in enumerate(ftp_urls):
                ftp_url = ftp_url.strip()
                if not ftp_url:
                    continue
                md5 = md5s[i].strip() if i < len(md5s) and md5s[i].strip() else None
                sample_to_files.setdefault(sample_acc, []).append((ftp_url, md5))

        return sample_to_files

    # === Download helpers ===
    def download_from_ftp(self, ftp_url: str, destination: Path) -> Path:
        """
        FTP プロトコルを用いてファイルをダウンロードします。再試行機構付き。

        Parameters
        ----------
        ftp_url : str
            ファイルの FTP URL。スキームがない場合は自動的に付与されます。
        destination : Path
            保存先のファイルパス。

        Returns
        -------
        Path
            保存先パス
        """
        # スキームが無い場合は付与
        if not ftp_url.startswith("ftp://"):
            ftp_url = "ftp://" + ftp_url
        return self._retry(self._download_ftp_once, ftp_url, destination)

    def _download_ftp_once(self, ftp_url: str, destination: Path) -> Path:
        """FTP ダウンロードを 1 回試行する。"""
        parse = urlparse(ftp_url)
        ftp_server = parse.netloc
        ftp_path = parse.path
        filename = os.path.basename(ftp_path)
        logger.info(f"{filename} を {destination} にFTP経由でダウンロードします")
        # FTP コネクションを開いてファイルを取得
        with FTP(ftp_server) as ftp:
            ftp.login()
            ftp.cwd(os.path.dirname(ftp_path))
            with open(destination, "wb") as fh:
                ftp.retrbinary("RETR " + filename, fh.write)
        logger.info(f"{filename} を {destination} にダウンロードしました")
        return destination

    def download_from_http(self, ftp_url: str, destination: Path) -> Path:
        """
        HTTP 経由でファイルをダウンロードします。FTP URL を HTTP に変換します。

        Parameters
        ----------
        ftp_url : str
            元々の FTP URL (またはスキームのないパス)。先頭を http:// に変換します。
        destination : Path
            保存先のファイルパス

        Returns
        -------
        Path
            保存先パス
        """
        # ftp:// → http:// へ変換
        if ftp_url.startswith("ftp://"):
            http_url = "http://" + ftp_url[len("ftp://") :]
        else:
            http_url = "http://" + ftp_url
        return self._retry(self._download_http_once, http_url, destination)

    def _download_http_once(self, http_url: str, destination: Path) -> Path:
        """HTTP ダウンロードを 1 回試行する。"""
        filename = os.path.basename(urlparse(http_url).path)
        logger.info(f"{filename} を {destination} にHTTP経由でダウンロードします")
        # ストリーミングでダウンロードしながら保存
        with requests.get(http_url, stream=True, timeout=20) as resp:
            resp.raise_for_status()
            with open(destination, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fh.write(chunk)
        logger.info(f"{filename} を {destination} にダウンロードしました")
        return destination

    def _retry(self, func, url: str, destination: Path) -> Path:
        """
        汎用の再試行ラッパー。指定回数失敗すると例外を送出します。

        Parameters
        ----------
        func : callable
            1 回のダウンロードを行う関数。
        url : str
            ダウンロード対象の URL。
        destination : Path
            保存先パス。

        Returns
        -------
        Path
            ダウンロードに成功した場合は destination を返します。
        """
        for attempt in range(1, self.max_retries + 1):
            try:
                return func(url, destination)
            except Exception as exc:
                if attempt >= self.max_retries:
                    logger.error(f"{url} のダウンロードに失敗しました: {exc}")
                    raise
                # 次の試行まで待機
                wait_seconds = 5
                logger.warning(
                    f"{url} のダウンロードに失敗しました (試行 {attempt}/{self.max_retries}): {exc}。"
                    f"{wait_seconds} 秒後に再試行します。"
                )
                time.sleep(wait_seconds)

    def download_sample_data(self, sample_acc: str, ftp_urls: list) -> list:
        """
        1 つのサンプルに対し関連する FASTQ ファイルをすべてダウンロードします。

        Parameters
        ----------
        sample_acc : str
            サンプルアクセッション。
        ftp_urls : list[str]
            そのサンプルに紐づく FTP URL のリスト。

        Returns
        -------
        list[Path]
            ダウンロードした FASTQ ファイルのパス一覧。
        """
        sample_dir = self.config.raw_data_dir / sample_acc
        sample_dir.mkdir(exist_ok=True)
        fastq_files: list = []
        # プロトコルに応じて使用するダウンロード関数を選択
        if self.protocol == "ftp":
            download_func = self.download_from_ftp
        else:
            download_func = self.download_from_http

        # ThreadPoolExecutor を用いた並列ダウンロード
        with ThreadPoolExecutor(max_workers=self.config.args.workers) as executor:
            future_to_url: dict = {}
            for url in ftp_urls:
                # FASTQ/GZ ファイルのみ対象
                if not url.endswith(".gz"):
                    continue
                # 保存先ファイル名を決定
                full_ftp_url = url if url.startswith("ftp://") else f"ftp://{url}"
                path = urlparse(full_ftp_url).path
                dest_path = sample_dir / os.path.basename(path)
                # ダウンロードタスクを送信
                future = executor.submit(download_func, url, dest_path)
                future_to_url[future] = url

            for future in as_completed(future_to_url):
                url = future_to_url[future]
                try:
                    dest = future.result()
                    fastq_files.append(dest)
                except Exception as exc:
                    # 既に _retry 内でエラーログを出力しているが、ここでも通知
                    logger.error(f"{url} からダウンロードに失敗しました: {exc}")
        return fastq_files