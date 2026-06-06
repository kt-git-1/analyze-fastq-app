import argparse
import logging
import shutil
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Optional
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

from config import TqdmLoggingHandler
from modules.ena_downloader import ENADownloader, verify_file_md5


def shorten_middle(text: str, width: int) -> str:
    """長い表示名を中央省略して、進捗バーの折り返しを避ける。"""
    if width <= 0:
        return ""
    if len(text) <= width:
        return text
    if width <= 3:
        return "." * width
    left = (width - 3) // 2
    right = width - 3 - left
    return text[:left] + "..." + text[-right:]


def download_display_name(sample_acc: str, filename: str, terminal_width: Optional[int] = None) -> str:
    """現在ファイル用の短い表示名を作る。"""
    if terminal_width is None:
        terminal_width = shutil.get_terminal_size((100, 20)).columns
    # tqdm の速度・サイズ表示が右側に出るため、desc は控えめに固定する。
    max_width = max(24, min(56, terminal_width // 2))
    sample_prefix = shorten_middle(sample_acc, 18)
    filename_width = max(12, max_width - len(sample_prefix) - 3)
    return "%s / %s" % (sample_prefix, shorten_middle(filename, filename_width))


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
    expected_md5: Optional[str] = None,
    progress_enabled: bool = True,
    progress_position: int = 1,
    display_name: Optional[str] = None,
) -> Path:
    """HTTPS ダウンロード（再開対応・MD5 検証付き）"""
    destination.parent.mkdir(parents=True, exist_ok=True)
    visible_name = display_name or shorten_middle(destination.name, 48)

    headers = {}
    mode = "wb"

    resume = False
    if destination.exists() and expected_md5 is None:
        size = destination.stat().st_size
        if size > 0:
            headers["Range"] = f"bytes={size}-"
            mode = "ab"
            resume = True
            logging.info("ダウンロードを再開します: %s (%d bytes 取得済み)", visible_name, size)
    elif destination.exists() and expected_md5 is not None:
        if verify_file_md5(destination, expected_md5):
            logging.info("MD5が一致したためスキップします: %s", visible_name)
            return destination
        logging.warning("MD5が一致しないため最初から再ダウンロードします: %s", visible_name)
        logging.debug("MD5不一致ファイル: %s", destination.name)
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
                    logging.warning("サーバーが再開ダウンロードに対応していないため、最初から再ダウンロードします: %s", visible_name)
                    destination.unlink(missing_ok=True)
                    headers.pop("Range", None)
                    mode = "wb"

                r.raise_for_status()

                total_size = int(r.headers.get("Content-Length", 0))
                initial = destination.stat().st_size if mode == "ab" and destination.exists() else 0

                with open(destination, mode) as f, tqdm(
                    total=total_size + initial if total_size else None,
                    initial=initial,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=visible_name,
                    position=progress_position,
                    leave=False,
                    dynamic_ncols=True,
                    disable=not progress_enabled,
                ) as pbar:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))

            if expected_md5 and not verify_file_md5(destination, expected_md5):
                if attempt < 5:
                    logging.warning(
                        "ダウンロード後のMD5が一致しません (%d/5): %s。再ダウンロードします",
                        attempt,
                        visible_name,
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
                        visible_name,
                    )
                    logging.debug("MD5不一致を保持したファイル: %s", destination.name)

            logging.info("ダウンロードが完了しました: %s", visible_name)
            return destination

        except Exception as e:
            last_err = e
            wait = min(60, 2**attempt)
            logging.warning(
                "ダウンロードに失敗しました (%d/5): %s。%d秒待って再試行します。詳細: %s",
                attempt,
                visible_name,
                wait,
                e,
            )
            time.sleep(wait)

    raise RuntimeError(f"ダウンロードに失敗しました: {url}") from last_err


def main():
    parser = argparse.ArgumentParser(description="ENA downloader (HTTPS CLI for testing)")
    parser.add_argument("--project", required=True, help="ENA project accession (e.g. PRJEB19970)")
    parser.add_argument("--out", required=True, help="base output dir (e.g. /path/to/raw_data)")
    parser.add_argument("--workers", type=int, default=1, help="download workers (default: 1)")
    parser.add_argument("--no-progress", action="store_true", help="disable terminal progress bars")
    args = parser.parse_args()
    progress_enabled = not args.no_progress

    handler = TqdmLoggingHandler() if progress_enabled else logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    logging.basicConfig(
        level=logging.INFO,
        handlers=[handler],
    )

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

    logging.info("検出されたサンプル数: %d", len(sample_to_urls))

    session_dl = build_https_session()

    total_files = sum(
        1 for pairs in sample_to_urls.values() for u, _ in pairs if u.endswith(".gz")
    )
    file_idx = 0
    if progress_enabled:
        logging.info(
            "ENA download: %s | サンプル: %d | ファイル: %d | ダウンロード並列数: %d",
            args.project,
            len(sample_to_urls),
            total_files,
            args.workers,
        )

    with tqdm(
        total=total_files,
        desc="全体進捗",
        unit="file",
        position=0,
        dynamic_ncols=True,
        disable=not progress_enabled,
    ) as overall_pbar:
        for sample_acc, url_md5_pairs in sample_to_urls.items():
            sample_dir = project_dir / sample_acc
            sample_dir.mkdir(parents=True, exist_ok=True)

            done_flag = sample_dir / ".done"
            sample_file_count = sum(1 for u, _ in url_md5_pairs if u.endswith(".gz"))
            if done_flag.exists():
                logging.info("スキップ (完了済み): %s", sample_acc)
                file_idx += sample_file_count
                overall_pbar.update(sample_file_count)
                continue

            for u, md5 in url_md5_pairs:
                if not u.endswith(".gz"):
                    continue
                file_idx += 1

                https_url = to_https_url(u)
                filename = Path(urlparse(https_url).path).name
                dest = sample_dir / filename

                short_current = download_display_name(sample_acc, filename)
                overall_pbar.set_postfix_str(
                    "%d/%d files | %s" % (file_idx, total_files, shorten_middle(sample_acc, 18)),
                    refresh=False,
                )
                logging.info(
                    "ダウンロード開始: %d/%d | %s",
                    file_idx,
                    total_files,
                    short_current,
                )
                logging.debug("ダウンロード対象: %s / %s", sample_acc, filename)
                download_via_https(
                    session_dl,
                    https_url,
                    dest,
                    expected_md5=md5,
                    progress_enabled=progress_enabled,
                    progress_position=1,
                    display_name=short_current,
                )
                overall_pbar.update(1)

            done_flag.touch()
            logging.info("完了フラグを記録しました: %s", sample_acc)


if __name__ == "__main__":
    main()
