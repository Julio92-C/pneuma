import csv
import pytest
from conftest import DATA, load_bin

ar = load_bin("arrange_reports")


def lines(name):
    return (DATA / name).read_text().splitlines()


def test_read_from_perc_strips_preamble():
    out = ar.read_from_perc(lines("kraken_combined.txt"))
    assert out[0].startswith("#perc\t")
    assert len(out) == 5


def test_read_from_perc_requires_one_header():
    with pytest.raises(ValueError):
        ar.read_from_perc(["#Number of Samples: 1", "no header here"])


def test_arrange_bracken_shape_and_values():
    rows = ar.arrange_bracken(ar.read_from_perc(lines("bracken_combined.txt")))
    assert rows[0] == ["#perc", "tot_all", "S1", "S2", "taxid", "name"]
    by_name = {r[-1]: r for r in rows[1:]}
    assert by_name["Streptococcus pneumoniae"][1:4] == ["60", "30", "30"]
    # denominator = root(90) + unclassified(10) = 100
    assert by_name["Streptococcus pneumoniae"][0] == "60"
    assert by_name["Haemophilus influenzae"][0] == "30"
    assert by_name["Haemophilus influenzae"][4] == "727"


def test_cli_writes_both_files(tmp_path):
    k_out, b_out = tmp_path / "k.txt", tmp_path / "b.csv"
    ar.main(["--kraken", str(DATA / "kraken_combined.txt"),
             "--bracken", str(DATA / "bracken_combined.txt"),
             "--kraken-out", str(k_out), "--bracken-out", str(b_out)])
    assert k_out.read_text().startswith("#perc\t")
    with open(b_out) as fh:
        assert next(csv.reader(fh))[0] == "#perc"


REAL_STYLE = [
    "#Number of Samples: 2",
    "#Total Number of Reads: 200",
    "#EBF0-0001\treports/a.txt",
    "#SPF0-0001\treports/b.txt",
    "#perc\ttot_all\ttot_lvl\t1_all\t1_lvl\t2_all\t2_lvl\tlvl_type\ttaxid\tname",
    "0.0000\t0\t0\tU\t0\tunclassified",
    "100.0000\t200\t0\t100\t0\t100\t0\tR\t1\troot",
    "100.0000\t200\t200\t100\t100\t100\t100\tS\t2697049\t  Severe acute respiratory syndrome coronavirus 2",
]


def test_real_combine_kreports_style():
    """Numbered columns + totals-only unclassified row, as written by KrakenTools 1.2."""
    names = ar.sample_names_from_preamble(REAL_STYLE)
    assert names == ["EBF0-0001", "SPF0-0001"]
    rows = ar.arrange_bracken(ar.read_from_perc(REAL_STYLE), names)
    assert rows[0] == ["#perc", "tot_all", "EBF0-0001", "SPF0-0001", "taxid", "name"]
    assert rows[1][:4] == ["0", "0", "0", "0"]
    assert rows[3][1:4] == ["200", "100", "100"]
