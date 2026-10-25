"""mixed_strain_check.py: minor-allele sites per kb on pooled reads vs the pooled assembly."""
import csv
import gzip

from conftest import DATA, load_bin

mx = load_bin("mixed_strain_check")
C = DATA / "pathogen_contracts"


def run(tmp_path, pileup, **kw):
    out = tmp_path / "u.mixed.tsv"
    args = {"min-depth": "8", "min-maf": "0.25", "min-alt": "3", "max-het-per-kb": "1.0", "min-evaluable": "1000"}
    args.update({k.replace("_", "-"): str(v) for k, v in kw.items()})
    argv = ["--unit", "u", "--mpileup", str(pileup), "--out", str(out)]
    for k, v in args.items():
        argv += ["--" + k, v]
    assert mx.main(argv) == 0
    return list(csv.DictReader(open(out), delimiter="\t"))[0]


def test_single_strain(tmp_path):
    row = run(tmp_path, C / "mpileup_single.txt")
    assert (row["evaluable_sites"], row["het_sites"], row["verdict"]) == ("2000", "0", "single")


def test_two_strains(tmp_path):
    row = run(tmp_path, C / "mpileup_mixed.txt")
    assert row["het_sites"] == "100" and row["het_per_kb"] == "50.000" and row["verdict"] == "mixed"
    assert list(row) == mx.OUT_COLUMNS


def test_gzip_input_and_not_evaluable(tmp_path):
    gz = tmp_path / "p.mpileup.gz"
    gz.write_bytes(gzip.compress((C / "mpileup_single.txt").read_bytes()))
    assert run(tmp_path, gz, min_evaluable=5000)["verdict"] == "not_evaluable"


def test_base_string_parsing():
    counts = mx.count_bases("A", ".,^~.$+2CG.-1T,*Gg>")
    assert counts["A"] == [3, 2]          # '.', '.', '.' forward; ',', ',' reverse (indel bases skipped)
    assert counts["G"] == [1, 1] and counts["C"] == [0, 0]


def test_single_strand_alt_is_not_het():
    lines = ["c\t1\tA\t10\t......GGGG\tIIIIIIIIII"]
    assert mx.check(lines, 8, 0.25, 3) == (1, 0)
