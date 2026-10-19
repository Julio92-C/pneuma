import csv
import re

import pytest

from conftest import load_bin

m = load_bin("rcf_for_meridian")


def _triplet(name):
    return [name, "", ""]


def make_raw(path):
    """Raw Recentrifuge 2.1.2 layout with anonymised IDs (ragged taxid rows)."""
    samples = ["S1", "S2", "CTRL1"]
    h = ["Samples"]
    st = ["Stats"]
    for s in samples:
        h += _triplet(f"kraken/{s}.kraken2.classifiedreads")
        st += ["cnt", "una", "sco"]
    for s in samples:
        h += _triplet(f"kraken/{s}.kraken2.classifiedreads_EXCLUSIVE_species")
        st += ["cnt", "una", "sco"]
    h += _triplet("kraken/S1.kraken2.classifiedreads_SHARED_genus")
    st += ["cnt", "una", "sco"]
    h += _triplet("SHARED_SUMMARY") + ["Details", "Details", "Name", "Rank"]
    st += ["cnt", "una", "sco", "Details", "Details", "Name", "Rank"]
    idrow = ["Id"] + [""] * (len(h) - 1)
    r1 = ["562"] + [str(i) for i in range(1, len(h) - 2)] + ["Escherichia coli", "species"]
    r2 = ["9606", "1", "2", "3"]  # ragged
    with open(path, "w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows([h, st, idrow, r1, r2])


def test_trim_keeps_one_triplet_per_sample(tmp_path):
    raw, out = tmp_path / "raw.csv", tmp_path / "out.csv"
    make_raw(raw)
    m.main(["--rcf", str(raw), "--samples", "S1,S2,CTRL1", "--out", str(out)])
    rows = list(csv.reader(open(out)))
    assert all(len(r) == 10 for r in rows)
    assert rows[0] == ["Samples", "kraken/S1", "", "", "kraken/S2", "", "", "kraken/CTRL1", "", ""]
    assert rows[1] == ["Stats"] + ["cnt", "una", "sco"] * 3
    assert rows[2][0] == "Id"
    assert rows[3][:4] == ["562", "1", "2", "3"]
    assert rows[4] == ["9606", "1", "2", "3"] + [""] * 6  # ragged row padded
    text = open(out).read()
    assert "EXCLUSIVE" not in text and "SHARED" not in text and "Details" not in text
    # MERIDIAN name stripping yields unique, exact sample IDs
    stripped = [re.sub(r"\..*", "", re.sub(r"^.*?/", "", c)) for c in rows[0][1:] if c]
    assert stripped == ["S1", "S2", "CTRL1"] and len(set(stripped)) == 3


def test_order_follows_samples_argument(tmp_path):
    raw, out = tmp_path / "raw.csv", tmp_path / "out.csv"
    make_raw(raw)
    m.main(["--rcf", str(raw), "--samples", "CTRL1,S1", "--out", str(out)])
    assert list(csv.reader(open(out)))[0][1::3] == ["kraken/CTRL1", "kraken/S1"]


def test_missing_sample_errors(tmp_path):
    raw, out = tmp_path / "raw.csv", tmp_path / "out.csv"
    make_raw(raw)
    with pytest.raises(SystemExit, match="'S9'"):
        m.main(["--rcf", str(raw), "--samples", "S1,S9", "--out", str(out)])


def test_duplicate_triplet_errors(tmp_path):
    raw, out = tmp_path / "raw.csv", tmp_path / "out.csv"
    rows = [["Samples"] + ["kraken/S1.kraken2.classifiedreads", "", ""] * 2,
            ["Stats"] + ["cnt", "una", "sco"] * 2, ["Id"] + [""] * 6]
    with open(raw, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    with pytest.raises(SystemExit, match="exactly 1"):
        m.main(["--rcf", str(raw), "--samples", "S1", "--out", str(out)])
