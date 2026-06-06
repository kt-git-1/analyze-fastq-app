import argparse
import logging
import shutil
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Callable, List, Optional
from urllib.parse import urlparse

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from modules.ena_downloader import ENADownloader, verify_file_md5


ProgressCallback = Callable[[int, Optional[int]], None]


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


def format_bytes(size: Optional[int]) -> str:
    if size is None:
        return "不明"
    value = float(size)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            if unit == "B":
                return "%d%s" % (int(value), unit)
            return "%.1f%s" % (value, unit)
        value /= 1024
    return "%.1fTB" % value


def format_rate(bytes_per_second: float) -> str:
    if bytes_per_second <= 0:
        return "-"
    return "%s/s" % format_bytes(int(bytes_per_second))


def format_eta(seconds: Optional[float]) -> str:
    if seconds is None:
        return "計算中"
    seconds = max(0, int(seconds))
    minutes, sec = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return "約 %d時間%d分" % (hours, minutes)
    if minutes:
        return "約 %d分%d秒" % (minutes, sec)
    return "約 %d秒" % sec


def progress_bar(done: int, total: int, width: int = 30) -> str:
    if total <= 0:
        return "[" + "-" * width + "]"
    filled = int(width * min(done, total) / total)
    return "[" + "#" * filled + "-" * (width - filled) + "]"


class DownloadDashboard:
    """ENA download用の固定表示ダッシュボード。"""

    def __init__(self, project: str, total_samples: int, total_files: int, parallel: int, enabled: bool = True):
        self.project = project
        self.total_samples = total_samples
        self.total_files = total_files
        self.parallel = parallel
        self.enabled = enabled
        self.file_index = 0
        self.sample_acc = "-"
        self.sample_done = 0
        self.sample_total = 0
        self.filename = "-"
        self.downloaded_bytes = 0
        self.total_bytes = None
        self.file_started_at = None
        self.events: List[str] = []
        self._rendered_lines = 0
        self._last_rendered_at = 0.0

    def add_event(self, message: str) -> None:
        timestamp = time.strftime("%H:%M:%S")
        self.events.append("%s  %s" % (timestamp, message))
        self.events = self.events[-3:]
        self.render(force=True)

    def start_file(
        self,
        file_index: int,
        sample_acc: str,
        sample_done: int,
        sample_total: int,
        filename: str,
    ) -> None:
        self.file_index = file_index
        self.sample_acc = sample_acc
        self.sample_done = sample_done
        self.sample_total = sample_total
        self.filename = filename
        self.downloaded_bytes = 0
        self.total_bytes = None
        self.file_started_at = time.monotonic()
        self.add_event("ダウンロード開始: %s / %s" % (sample_acc, shorten_middle(filename, 42)))

    def update_file(self, downloaded_bytes: int, total_bytes: Optional[int]) -> None:
        self.downloaded_bytes = downloaded_bytes
        self.total_bytes = total_bytes
        self.render()

    def finish_file(self, md5_checked: bool = False) -> None:
        self.sample_done += 1
        if md5_checked:
            self.add_event("MD5確認OK: %s" % shorten_middle(self.filename, 42))
        else:
            self.add_event("ダウンロード完了: %s" % shorten_middle(self.filename, 42))
        self.add_event("次のファイルへ移動")

    def skip_sample(self, sample_acc: str, sample_file_count: int) -> None:
        self.file_index += sample_file_count
        self.sample_acc = sample_acc
        self.sample_done = sample_file_count
        self.sample_total = sample_file_count
        self.add_event("スキップ (完了済み): %s" % sample_acc)

    def render_text(self) -> str:
        percent = int(round(self.file_index / self.total_files * 100)) if self.total_files else 100
        now = time.monotonic()
        elapsed = 0.0 if self.file_started_at is None else max(0.0, now - self.file_started_at)
        rate = self.downloaded_bytes / elapsed if elapsed > 0 else 0.0
        eta = None
        if self.total_bytes and rate > 0:
            eta = (self.total_bytes - self.downloaded_bytes) / rate

        terminal_width = shutil.get_terminal_size((100, 20)).columns
        filename_width = max(24, min(54, terminal_width - 12))
        sample_width = max(18, min(36, terminal_width - 32))
        event_lines = self.events or ["-"]

        return "\n".join(
            [
                "ENA download: %s" % self.project,
                "",
                "全体進捗  %s  %3d%%  %d/%d files" % (
                    progress_bar(self.file_index, self.total_files),
                    percent,
                    self.file_index,
                    self.total_files,
                ),
                "サンプル  %-*s  %d/%d files" % (
                    sample_width,
                    shorten_middle(self.sample_acc, sample_width),
                    self.sample_done,
                    self.sample_total,
                ),
                "現在      %s" % shorten_middle(self.filename, filename_width),
                "サイズ    %s / %s" % (format_bytes(self.downloaded_bytes), format_bytes(self.total_bytes)),
                "速度      %s" % format_rate(rate),
                "残り時間  %s" % format_eta(eta),
                "",
                "最近のイベント",
            ]
            + ["  %s" % line for line in event_lines]
        )

    def render(self, force: bool = False) -> None:
        if not self.enabled:
            return
        now = time.monotonic()
        if not force and now - self._last_rendered_at < 0.2:
            return
        text = self.render_text()
        if self._rendered_lines:
            sys.stderr.write("\033[%dF\033[J" % self._rendered_lines)
        sys.stderr.write(text + "\n")
        sys.stderr.flush()
        self._rendered_lines = text.count("\n") + 1
        self._last_rendered_at = now

    def close(self) -> None:
        self.render(force=True)
        if self.enabled:
            sys.stderr.write("\n")
            sys.stderr.flush()


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
    display_name: Optional[str] = None,
    progress_callback: Optional[ProgressCallback] = None,
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
                visible_total = total_size + initial if total_size else None
                downloaded = initial

                if progress_callback:
                    progress_callback(downloaded, visible_total)

                with open(destination, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            if progress_callback:
                                progress_callback(downloaded, visible_total)

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

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
    handler.setLevel(logging.WARNING if progress_enabled else logging.INFO)
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
    dashboard = DownloadDashboard(
        args.project,
        total_samples=len(sample_to_urls),
        total_files=total_files,
        parallel=args.workers,
        enabled=progress_enabled,
    )
    dashboard.render(force=True)
    try:
        for sample_acc, url_md5_pairs in sample_to_urls.items():
            sample_dir = project_dir / sample_acc
            sample_dir.mkdir(parents=True, exist_ok=True)

            done_flag = sample_dir / ".done"
            sample_file_count = sum(1 for u, _ in url_md5_pairs if u.endswith(".gz"))
            if done_flag.exists():
                logging.info("スキップ (完了済み): %s", sample_acc)
                dashboard.skip_sample(sample_acc, sample_file_count)
                file_idx += sample_file_count
                continue

            sample_done = 0
            for u, md5 in url_md5_pairs:
                if not u.endswith(".gz"):
                    continue
                file_idx += 1

                https_url = to_https_url(u)
                filename = Path(urlparse(https_url).path).name
                dest = sample_dir / filename

                short_current = download_display_name(sample_acc, filename)
                logging.info(
                    "ダウンロード開始: %d/%d | %s",
                    file_idx,
                    total_files,
                    short_current,
                )
                logging.debug("ダウンロード対象: %s / %s", sample_acc, filename)
                dashboard.start_file(file_idx, sample_acc, sample_done, sample_file_count, filename)
                download_via_https(
                    session_dl,
                    https_url,
                    dest,
                    expected_md5=md5,
                    progress_enabled=progress_enabled,
                    display_name=short_current,
                    progress_callback=dashboard.update_file if progress_enabled else None,
                )
                sample_done += 1
                dashboard.finish_file(md5_checked=bool(md5))

            done_flag.touch()
            logging.info("完了フラグを記録しました: %s", sample_acc)
    finally:
        dashboard.close()


if __name__ == "__main__":
    main()
