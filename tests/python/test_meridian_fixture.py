"""The synthetic MERIDIAN fixture must stay reproducible and PNEUMA-realistic."""
import csv
import importlib.util
import filecmp
from pathlib import Path

from conftest import DATA, load_bin

FIX = DATA / "meridian_fixture"
SAMPLES = ["CTRL1"] + [f"S0{i}" for i in range(1, 9)]


def _gen():
    spec = importlib.util.spec_from_file_location(
        "make_meridian_fixture", Path(__file__).resolve().parents[1] / "make_meridian_fixture.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_generator_is_deterministic_and_committed_output_current(tmp_path):
    _gen().main([str(tmp_path)])
    for rel in ["kraken2_combined_reports.txt", "bracken_arranged.csv", "pneuma.rcf.data.csv",
                "summary_reportClean.csv", "metadata.csv", "raw/kraken_combined_raw.txt",
                "raw/bracken_combined_raw.txt"]:
        assert filecmp.cmp(tmp_path / rel, FIX / rel, shallow=False), f"{rel} differs from generator output"


def test_arrange_reports_reproduces_arranged_files(tmp_path):
    ar = load_bin("arrange_reports")
    k, b = tmp_path / "k.txt", tmp_path / "b.csv"
    ar.main(["--kraken", str(FIX / "raw/kraken_combined_raw.txt"),
             "--bracken", str(FIX / "raw/bracken_combined_raw.txt"),
             "--kraken-out", str(k), "--bracken-out", str(b)])
    assert k.read_text() == (FIX / "kraken2_combined_reports.txt").read_text()
    assert b.read_text() == (FIX / "bracken_arranged.csv").read_text()


def test_rcf_trims_to_one_triplet_per_sample(tmp_path):
    rcf = load_bin("rcf_for_meridian")
    out = tmp_path / "o.csv"
    rcf.main(["--rcf", str(FIX / "pneuma.rcf.data.csv"), "--samples", ",".join(SAMPLES), "--out", str(out)])
    rows = list(csv.reader(open(out)))
    assert len(rows[0]) == 1 + 3 * len(SAMPLES)
    assert rows[0][0] == "Samples" and rows[1][0] == "Stats"


def test_abricate_layout_and_zero_hit_sample():
    merge = load_bin("merge_abricate")
    lines = (FIX / "summary_reportClean.csv").read_text().splitlines()
    assert lines[0] == merge.HEADER
    rows = [l.split(",") for l in lines[1:]]
    assert all(len(r) == 15 for r in rows)
    assert {r[0].split("_kraken2ID")[0] for r in rows} == set(SAMPLES) - {"S07"}
    assert len({r[5] for r in rows if r[11] == "vfdb"}) >= 8
    assert all("_kraken:taxid|" in r[1] for r in rows)


def test_metadata_has_groups_and_control():
    rows = list(csv.DictReader(open(FIX / "metadata.csv")))
    assert [r["sample"] for r in rows] == SAMPLES[1:] + ["CTRL1"]
    assert {r["group"] for r in rows if r["sample"] != "CTRL1"} == {"A", "B"}
