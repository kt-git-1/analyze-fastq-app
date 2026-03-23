import argparse
import logging
import time
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

from modules.ena_downloader import ENADownloader, verify_file_md5


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
    expected_md5: str | None = None,
) -> Path:
    """HTTPS ダウンロード（再開対応・MD5 検証付き）"""
    destination.parent.mkdir(parents=True, exist_ok=True)

    headers = {}
    mode = "wb"

    resume = False
    if destination.exists() and expected_md5 is None:
        size = destination.stat().st_size
        if size > 0:
            headers["Range"] = f"bytes={size}-"
            mode = "ab"
            resume = True
            logging.info(f"再開DL: {destination.name} ({size} bytes)")
    elif destination.exists() and expected_md5 is not None:
        if verify_file_md5(destination, expected_md5):
            logging.info(f"MD5一致 (スキップ): {destination.name}")
            return destination
        logging.warning(f"MD5不一致 → 最初から再DL: {destination.name}")
        destination.unlink(missing_ok=True)

    last_err = None
    for attempt in range(1, 6):
        try:
            with session.get(
                url,
                stream=True,
                headers=headers,
                timeout=(10, 300),
            ) as r:
                if "Range" in headers and r.status_code == 200:
                    logging.warning(f"Range無視 → 最初から再DL: {destination.name}")
                    destination.unlink(missing_ok=True)
                    headers.pop("Range", None)
                    mode = "wb"

                r.raise_for_status()

                total_size = int(r.headers.get("Content-Length", 0))
                initial = destination.stat().st_size if mode == "ab" and destination.exists() else 0

                with (
                    open(destination, mode) as f,
                    tqdm(
                        total=total_size + initial if total_size else None,
                        initial=initial,
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                        desc=destination.name,
                        leave=False,
                    ) as pbar,
                ):
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

            if expected_md5 and not verify_file_md5(destination, expected_md5):
                if attempt < 5:
                    logging.warning(
                        f"DL後 MD5不一致 ({attempt}/5): {destination.name} → 再ダウンロード"
                    )
                    destination.unlink(missing_ok=True)
                    headers.pop("Range", None)
                    mode = "wb"
                    resume = False
                    continue
                else:
                    logging.warning(
                        "MD5不一致が解消しません: %s "
                        "(ENA メタデータと実ファイルの不整合の可能性あり → ファイルを保持して続行)",
                        destination.name,
                    )

            logging.info(f"DL完了: {destination.name}")
            return destination

        except Exception as e:
            last_err = e
            wait = min(60, 2**attempt)
            logging.warning(
                f"DL失敗({attempt}/5): {destination.name} : {e} → {wait}s待機"
            )
            time.sleep(wait)

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

    # ENA APIは本体を利用（返り値は {sample_acc: [(url, md5), ...]}）
    tsv = downloader.get_api_response(args.project, session_api)
    sample_to_urls = downloader.parse_response_with_checksums(tsv)

    logging.info(f"検出されたサンプル数: {len(sample_to_urls)}")

    session_dl = build_https_session()

    total_files = sum(
        1 for pairs in sample_to_urls.values() for u, _ in pairs if u.endswith(".gz")
    )
    file_idx = 0

    for sample_acc, url_md5_pairs in sample_to_urls.items():
        sample_dir = project_dir / sample_acc
        sample_dir.mkdir(parents=True, exist_ok=True)

        done_flag = sample_dir / ".done"
        if done_flag.exists():
            logging.info(f"SKIP (done): {sample_acc}")
            file_idx += sum(1 for u, _ in url_md5_pairs if u.endswith(".gz"))
            continue

        for u, md5 in url_md5_pairs:
            if not u.endswith(".gz"):
                continue
            file_idx += 1

            https_url = to_https_url(u)
            filename = Path(urlparse(https_url).path).name
            dest = sample_dir / filename

            logging.info(
                "[%d/%d] %s / %s", file_idx, total_files, sample_acc, filename,
            )
            download_via_https(session_dl, https_url, dest, expected_md5=md5)

        done_flag.touch()
        logging.info(f"MARK DONE: {sample_acc}")


if __name__ == "__main__":
    main()
