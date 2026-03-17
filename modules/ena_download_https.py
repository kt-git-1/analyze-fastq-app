import argparse
import logging
import time
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# 本体は変更しない（友人コード）
from modules.ena_downloader import ENADownloader


def build_https_session() -> requests.Session:
    """HTTPS + retry 用 Session"""
    session = requests.Session()

    retry = Retry(
        total=10,
        connect=10,
        read=10,
        status=10,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["GET", "HEAD"]),
        raise_on_status=False,
    )

    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=20,
        pool_maxsize=20,
    )
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def to_https_url(u: str) -> str:
    """ENAの ftp URL を https に変換"""
    u = u.strip()
    if u.startswith("ftp://"):
        return "https://" + u[len("ftp://") :]
    if u.startswith(("http://", "https://")):
        return u
    return "https://" + u.lstrip("/")


def download_via_https(
    session: requests.Session,
    url: str,
    destination: Path,
    chunk_size: int = 1024 * 1024,
) -> Path:
    """HTTPS ダウンロード（再開対応）"""
    destination.parent.mkdir(parents=True, exist_ok=True)

    headers = {}
    mode = "wb"

    if destination.exists():
        size = destination.stat().st_size
        if size > 0:
            headers["Range"] = f"bytes={size}-"
            mode = "ab"
            logging.info(f"再開DL: {destination.name} ({size} bytes)")

    last_err = None
    for attempt in range(1, 6):
        try:
            with session.get(
                url,
                stream=True,
                headers=headers,
                timeout=(10, 300),
            ) as r:
                # Range を無視されたら最初から
                if "Range" in headers and r.status_code == 200:
                    logging.warning(f"Range無視 → 最初から再DL: {destination.name}")
                    destination.unlink(missing_ok=True)
                    headers.pop("Range", None)
                    mode = "wb"

                r.raise_for_status()

                with open(destination, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)

            logging.info(f"DL完了: {destination.name}")
            return destination

        except Exception as e:
            last_err = e
            wait = min(60, 2**attempt)
            logging.warning(
                f"DL失敗({attempt}/5): {destination.name} : {e} → {wait}s待機"
            )
            time.sleep(wait)

    # last_err が None の可能性もあるので安全に例外化
    raise RuntimeError(f"download failed: {url}") from last_err


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    parser = argparse.ArgumentParser(description="ENA downloader (HTTPS CLI for testing)")
    parser.add_argument("--project", required=True, help="ENA project accession (e.g. PRJEB19970)")
    parser.add_argument("--out", required=True, help="base output dir (e.g. /path/to/raw_data)")
    parser.add_argument("--workers", type=int, default=1, help="download workers (default: 1)")
    args = parser.parse_args()

    # プロジェクト番号ディレクトリを作成
    project_dir = Path(args.out) / args.project
    project_dir.mkdir(parents=True, exist_ok=True)

    # 本体が期待する config 形式に合わせる（本体は変更しない）
    config = SimpleNamespace(
        raw_data_dir=project_dir,
        args=SimpleNamespace(workers=args.workers),
    )

    downloader = ENADownloader(config)
    session_api = build_https_session()

    # ENA APIは本体を利用（返り値は {sample_acc: [url, ...]}）
    tsv = downloader.get_api_response(args.project, session_api)
    sample_to_urls = downloader.parse_response_data(tsv)

    logging.info(f"検出されたサンプル数: {len(sample_to_urls)}")

    session_dl = build_https_session()

    # サンプル単位で保存: out/PROJECT/SAMPLE/ （main.py の --fastq_dir と互換）
    for sample_acc, urls in sample_to_urls.items():
        sample_dir = project_dir / sample_acc
        sample_dir.mkdir(parents=True, exist_ok=True)

        done_flag = sample_dir / ".done"
        if done_flag.exists():
            logging.info(f"SKIP (done): {sample_acc}")
            continue

        for u in urls:
            if not u.endswith(".gz"):
                continue

            https_url = to_https_url(u)
            filename = Path(urlparse(https_url).path).name
            dest = sample_dir / filename

            download_via_https(session_dl, https_url, dest)

        done_flag.touch()
        logging.info(f"MARK DONE: {sample_acc}")


if __name__ == "__main__":
    main()
