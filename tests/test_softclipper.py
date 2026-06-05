from types import SimpleNamespace

from modules.softclipper import SoftClipper


class FakeRead:
    def __init__(self, cigartuples, sequence="A" * 10, unmapped=False):
        self.cigartuples = cigartuples
        self.query_sequence = sequence
        self.query_qualities = [30] * len(sequence)
        self.is_unmapped = unmapped
        self.query_name = "read1"


def test_calculate_query_length_counts_query_consuming_ops(make_config):
    clipper = SoftClipper(make_config())

    assert clipper.calculate_query_length([(0, 10), (1, 2), (2, 3), (4, 5)]) == 17
    assert clipper.calculate_query_length(None) == 0


def test_process_read_leaves_unmapped_or_cigarless_reads_unchanged(make_config):
    clipper = SoftClipper(make_config())
    unmapped = FakeRead([(0, 10)], unmapped=True)
    cigarless = FakeRead(None)

    assert clipper.process_read(unmapped) is unmapped
    assert clipper.process_read(cigarless) is cigarless


def test_process_read_adds_five_base_soft_clip(make_config):
    read = FakeRead([(0, 10)], "A" * 10)

    out = SoftClipper(make_config()).process_read(read)

    assert out.cigartuples == [(4, 5), (0, 5)]
    assert out.query_sequence == "A" * 10


def test_process_read_merges_existing_soft_clip(make_config):
    read = FakeRead([(4, 2), (0, 10)], "A" * 12)

    out = SoftClipper(make_config()).process_read(read)

    assert out.cigartuples == [(4, 7), (0, 5)]


def test_process_read_preserves_deletions(make_config):
    read = FakeRead([(0, 8), (2, 3), (0, 2)], "A" * 10)

    out = SoftClipper(make_config()).process_read(read)

    assert out.cigartuples == [(4, 5), (0, 3), (2, 3), (0, 2)]


def test_process_read_returns_none_when_query_is_too_short(make_config):
    read = FakeRead([(0, 10)], "A" * 4)

    assert SoftClipper(make_config()).process_read(read) is None


def test_count_reads_returns_none_when_samtools_fails(monkeypatch, make_config, tmp_path):
    monkeypatch.setattr("modules.softclipper.subprocess.run", lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("bad")))

    assert SoftClipper(make_config())._count_reads(tmp_path / "reads.bam") is None
