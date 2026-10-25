"""coverage_gate.py: genome-wide depth and breadth from samtools coverage."""
import csv

import pytest

from conftest import DATA, load_bin

cg = load_bin("coverage_gate")
SNIPPET = DATA / "pathogen_contracts" / "samtools_coverage.txt"


def read(path):
    return list(csv.DictReader(open(path), delimiter="\t"))[0]


def test_captured_real_output(tmp_path):
    out = tmp_path / "u.gate.tsv"
    assert cg.main(["--unit", "MINI-S1.Streptococcus_pneumoniae", "--coverage", str(SNIPPET),
                    "--min-depth", "5", "--min-breadth", "0.5", "--out", str(out)]) == 0
    row = read(out)
    assert list(row) == cg.OUT_COLUMNS
    # captured line: MINI-TIGR4-SEG 1 120000 numreads 188, covbases 120000 (breadth 1.0), meandepth 8.61925
    assert row["mean_depth"] == "8.619" and row["breadth"] == "1.0000" and row["reads_mapped"] == "188"
    assert row["pass"] == "true"


def test_length_weighted_depth_and_breadth_over_chromosome_and_plasmid(tmp_path):
    cov = tmp_path / "c.txt"
    cov.write_text(SNIPPET.read_text().splitlines()[0] + "\n"
                   "chr\t1\t2000000\t3000\t1900000\t95\t12.5\t20.1\t59.8\n"
                   "plasmid\t1\t100000\t40\t20000\t20\t1\t19.5\t60\n")
    out = tmp_path / "u.gate.tsv"
    assert cg.main(["--unit", "u", "--coverage", str(cov), "--min-depth", "10", "--min-breadth", "0.5",
                    "--out", str(out)]) == 0
    row = read(out)
    # (12.5 * 2,000,000 + 1 * 100,000) / 2,100,000 = 11.952; (1,900,000 + 20,000) / 2,100,000 = 0.9143
    assert row["mean_depth"] == "11.952" and row["breadth"] == "0.9143" and row["reads_mapped"] == "3040"
    assert row["pass"] == "true"


@pytest.mark.parametrize("depth,breadth,expected", [
    ("8.7", "0.5", "false"),    # depth 8.619 < 8.7
    ("8.6", "0.5", "true"),
    ("8", "1.01", "false"),     # breadth 1.0 < 1.01 (the captured unit covers 100% of the reference)
])
def test_thresholds(tmp_path, depth, breadth, expected):
    out = tmp_path / "g.tsv"
    cg.main(["--unit", "u", "--coverage", str(SNIPPET), "--min-depth", depth, "--min-breadth", breadth,
             "--out", str(out)])
    assert read(out)["pass"] == expected


def test_no_reads_mapped_fails_the_gate(tmp_path):
    cov = tmp_path / "c.txt"
    cov.write_text(SNIPPET.read_text().splitlines()[0] + "\nchr\t1\t1000\t0\t0\t0\t0\t0\t0\n")
    out = tmp_path / "g.tsv"
    cg.main(["--unit", "u", "--coverage", str(cov), "--min-depth", "1", "--min-breadth", "0.1", "--out", str(out)])
    assert read(out)["pass"] == "false" and read(out)["mean_depth"] == "0.000"


def test_wrong_header_is_a_contract_error(tmp_path):
    cov = tmp_path / "c.txt"
    cov.write_text("rname\tstart\n")
    with pytest.raises(ValueError, match="samtools coverage"):
        cg.gate("u", str(cov), 1, 0.1)
