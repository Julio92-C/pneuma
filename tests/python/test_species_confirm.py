"""species_confirm.py: best skani hit must be the panel species with ANI >= ani_min."""
import csv

from conftest import DATA, load_bin

sc = load_bin("species_confirm")
SNIPPET = DATA / "pathogen_contracts" / "skani_dist.tsv"
REFS = DATA.parents[1] / "assets" / "panel_references.tsv"


def run(tmp_path, species, ani_min, skani=SNIPPET):
    out = tmp_path / "u.species.tsv"
    assert sc.main(["--unit", "u", "--species", species, "--ani-min", str(ani_min), "--skani", str(skani),
                    "--references", str(REFS), "--out", str(out)]) == 0
    return list(csv.DictReader(open(out), delimiter="\t"))[0]


def test_mitis_reads_labelled_pneumoniae_are_not_confirmed(tmp_path):
    row = run(tmp_path, "Streptococcus pneumoniae", 95)
    assert row["confirmed"] == "false" and row["reason"] == "best_hit_other_species"
    assert row["best_species"] == "Streptococcus mitis" and row["best_accession"] == "GCF_018603625.2"
    assert row["ani"] == "96.85"


def test_confirmed_and_ani_threshold(tmp_path):
    assert run(tmp_path, "Streptococcus mitis", 95)["confirmed"] == "true"
    low = run(tmp_path, "Streptococcus mitis", 97)
    assert low["confirmed"] == "false" and low["reason"] == "ani_below_min"


def test_empty_skani_output_is_not_confirmed(tmp_path):
    empty = tmp_path / "e.tsv"
    empty.write_text(SNIPPET.read_text().splitlines()[0] + "\n")
    row = run(tmp_path, "Streptococcus pneumoniae", 95, empty)
    assert row["confirmed"] == "false" and row["reason"] == "no_skani_hit" and row["ani"] == ""
    assert list(row) == sc.OUT_COLUMNS
