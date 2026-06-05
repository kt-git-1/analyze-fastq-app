import hashlib

import pytest

import modules.ena_download_https as https_module
from modules.ena_download_https import download_via_https, to_https_url


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
