import io
from conftest import DATA, load_bin

merge = load_bin("merge_abricate")


def test_single_header_and_rows_kept():
    buf = io.StringIO()
    n = merge.merge([DATA / "abricate_S1.csv", DATA / "abricate_S2_empty.csv"], buf)
    lines = buf.getvalue().splitlines()
    assert n == 2
    assert lines[0] == merge.HEADER
    assert sum(1 for l in lines if l.startswith("#FILE")) == 1
    assert lines[2].endswith("IncFII,")


def test_only_empty_inputs_gives_header_only():
    buf = io.StringIO()
    assert merge.merge([DATA / "abricate_S2_empty.csv"], buf) == 0
    assert buf.getvalue() == merge.HEADER + "\n"


def test_unexpected_header_raises(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("#FILE,SEQUENCE,GENE\nx,y,z\n")
    try:
        merge.merge([bad], io.StringIO())
    except ValueError as e:
        assert "header" in str(e)
    else:
        raise AssertionError("expected ValueError")
