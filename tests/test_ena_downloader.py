from types import SimpleNamespace

import pytest

from modules.ena_downloader import ENADownloader, verify_file_md5


def make_downloader(tmp_path, protocol="http", retries=2):
    raw_data_dir = tmp_path / "raw"
    raw_data_dir.mkdir(parents=True, exist_ok=True)
    config = SimpleNamespace(
        raw_data_dir=raw_data_dir,
        args=SimpleNamespace(workers=2, download_protocol=protocol, max_retries=retries),
    )
    return ENADownloader(config)


def test_verify_file_md5_matches_and_mismatches(tmp_path):
    path = tmp_path / "file.txt"
    path.write_text("abc")

    assert verify_file_md5(path, "900150983cd24fb0d6963f7d28e17f72") is True
    assert verify_file_md5(path, "bad") is False


def test_parse_response_with_checksums_handles_multiple_urls_and_missing_md5(tmp_path):
    data = (
        "sample_accession\tsubmitted_ftp\tsubmitted_md5\n"
        "S1\tftp.sra/a_1.fastq.gz;ftp.sra/a_2.fastq.gz\tmd51;md52\n"
        "S2\tftp.sra/b.fastq.gz\t\n"
    )

    parsed = make_downloader(tmp_path).parse_response_with_checksums(data)

    assert parsed == {
        "S1": [("ftp.sra/a_1.fastq.gz", "md51"), ("ftp.sra/a_2.fastq.gz", "md52")],
        "S2": [("ftp.sra/b.fastq.gz", None)],
    }


def test_parse_response_with_checksums_returns_empty_for_bad_header(tmp_path):
    assert make_downloader(tmp_path).parse_response_with_checksums("bad\theader\nx\ty") == {}


def test_download_from_http_converts_ftp_url(monkeypatch, tmp_path):
    downloader = make_downloader(tmp_path)
    calls = []
    monkeypatch.setattr(
        downloader,
        "_retry",
        lambda func, url, destination: calls.append((func, url, destination)) or destination,
    )

    dest = tmp_path / "out.fastq.gz"
    assert downloader.download_from_http("ftp://host/path/out.fastq.gz", dest) == dest
    assert calls[0][1] == "http://host/path/out.fastq.gz"


def test_retry_retries_then_succeeds(monkeypatch, tmp_path):
    downloader = make_downloader(tmp_path, retries=3)
    monkeypatch.setattr("modules.ena_downloader.time.sleep", lambda seconds: None)
    calls = {"count": 0}

    def flaky(url, destination):
        calls["count"] += 1
        if calls["count"] < 2:
            raise RuntimeError("temporary")
        return destination

    assert downloader._retry(flaky, "url", tmp_path / "dest") == tmp_path / "dest"
    assert calls["count"] == 2


def test_retry_raises_after_max_attempts(monkeypatch, tmp_path):
    downloader = make_downloader(tmp_path, retries=2)
    monkeypatch.setattr("modules.ena_downloader.time.sleep", lambda seconds: None)

    with pytest.raises(RuntimeError):
        downloader._retry(lambda url, destination: (_ for _ in ()).throw(RuntimeError("bad")), "url", tmp_path / "dest")


def test_download_sample_data_skips_non_gzip_and_uses_protocol(monkeypatch, tmp_path):
    downloader = make_downloader(tmp_path)
    calls = []

    def fake_download(url, destination):
        calls.append((url, destination))
        return destination

    monkeypatch.setattr(downloader, "download_from_http", fake_download)

    files = downloader.download_sample_data(
        "S1",
        ["ftp.sra/readme.txt", "ftp.sra/S1_R1.fastq.gz", "ftp://ftp.sra/S1_R2.fastq.gz"],
    )

    assert len(files) == 2
    assert {call[0] for call in calls} == {"ftp.sra/S1_R1.fastq.gz", "ftp://ftp.sra/S1_R2.fastq.gz"}
