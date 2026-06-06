import hashlib

import pytest

import modules.ena_download_https as https_module
from modules.ena_download_https import (
    DownloadDashboard,
    download_display_name,
    download_via_https,
    shorten_middle,
    to_https_url,
)


class FakeResponse:
    def __init__(self, chunks, status_code=206):
        self._chunks = chunks
        self.status_code = status_code
        self.headers = {"Content-Length": str(sum(len(c) for c in chunks))}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        return iter(self._chunks)


class FakeSession:
    def __init__(self, response):
        self.response = response
        self.calls = []

    def get(self, url, stream, headers, timeout):
        self.calls.append({"url": url, "headers": dict(headers)})
        return self.response


@pytest.mark.parametrize(
    ("url", "expected"),
    [
        ("ftp://host/path.fastq.gz", "https://host/path.fastq.gz"),
        ("http://host/path.fastq.gz", "http://host/path.fastq.gz"),
        ("https://host/path.fastq.gz", "https://host/path.fastq.gz"),
        ("host/path.fastq.gz", "https://host/path.fastq.gz"),
    ],
)
def test_to_https_url(url, expected):
    assert to_https_url(url) == expected


def test_shorten_middle_keeps_short_text_and_truncates_long_text():
    assert shorten_middle("short.fastq.gz", 20) == "short.fastq.gz"
    assert shorten_middle("abcdefghijklmnopqrstuvwxyz", 12) == "abcd...vwxyz"
    assert shorten_middle("abcdef", 3) == "..."


def test_download_display_name_limits_long_fastq_names():
    filename = "BER01_A__BER01_A_E16.1_user_TGTCTG__BER01_A_E16.1_user_TGTCTG_PC4T_20140922_hiseq3a.fastq.gz"

    display = download_display_name("SAMEA103910511", filename, terminal_width=90)

    assert display.startswith("SAMEA103910511 / ")
    assert "..." in display
    assert len(display) <= 45


def test_download_dashboard_renders_human_readable_status(monkeypatch):
    monkeypatch.setattr(https_module.time, "monotonic", lambda: 100.0)
    dashboard = DownloadDashboard("PRJEB19970", total_samples=14, total_files=216, parallel=4, enabled=False)
    dashboard.file_started_at = 90.0
    dashboard.start_file(
        153,
        "SAMEA103910511",
        7,
        12,
        "BER01_A__BER01_A_E16.1_user_TGTCTG_PC4T_20140922_hiseq3a.fastq.gz",
    )
    dashboard.update_file(420 * 1024 * 1024, 588 * 1024 * 1024)
    dashboard.finish_file(md5_checked=True)

    text = dashboard.render_text()

    assert "ENA download: PRJEB19970" in text
    assert "全体進捗" in text
    assert "153/216 files" in text
    assert "サンプル  SAMEA103910511" in text
    assert "8/12 files" in text
    assert "現在" in text
    assert "サイズ" in text
    assert "速度" in text
    assert "残り時間" in text
    assert "最近のイベント" in text
    assert "MD5確認OK" in text
    assert "root:" not in text
    assert "__main__:" not in text


def test_download_via_https_reports_progress(monkeypatch, tmp_path):
    dest = tmp_path / "file.fastq.gz"
    session = FakeSession(FakeResponse([b"ab", b"cd"], status_code=206))
    calls = []

    assert download_via_https(
        session,
        "https://example/file.fastq.gz",
        dest,
        progress_callback=lambda done, total: calls.append((done, total)),
    ) == dest

    assert calls == [(0, 4), (2, 4), (4, 4)]


def test_download_via_https_skips_when_existing_md5_matches(monkeypatch, tmp_path):
    dest = tmp_path / "file.fastq.gz"
    dest.write_bytes(b"complete")
    monkeypatch.setattr(https_module, "verify_file_md5", lambda path, expected: True)
    session = FakeSession(FakeResponse([b"new"]))

    assert download_via_https(session, "https://example/file.fastq.gz", dest, expected_md5="md5") == dest
    assert session.calls == []
    assert dest.read_bytes() == b"complete"


def test_download_via_https_resumes_and_restarts_when_range_is_ignored(tmp_path):
    dest = tmp_path / "file.fastq.gz"
    dest.write_bytes(b"partial")
    session = FakeSession(FakeResponse([b"fresh"], status_code=200))

    assert download_via_https(session, "https://example/file.fastq.gz", dest) == dest

    assert session.calls[0]["headers"] == {"Range": "bytes=7-"}
    assert dest.read_bytes() == b"fresh"


def test_download_via_https_redownloads_until_md5_matches(monkeypatch, tmp_path):
    dest = tmp_path / "file.fastq.gz"
    expected = hashlib.md5(b"ok").hexdigest()
    monkeypatch.setattr(https_module.time, "sleep", lambda seconds: None)
    session = FakeSession(FakeResponse([b"ok"]))

    assert download_via_https(session, "https://example/file.fastq.gz", dest, expected_md5=expected) == dest
    assert dest.read_bytes() == b"ok"
