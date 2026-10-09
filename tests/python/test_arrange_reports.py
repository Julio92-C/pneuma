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
