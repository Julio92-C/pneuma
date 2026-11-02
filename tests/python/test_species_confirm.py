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


def ecoli_shigella(tmp_path):
    p = tmp_path / "ecoli.skani.tsv"
    head = SNIPPET.read_text().splitlines()[0]
    p.write_text(head + "\n"
                 "GCF_002950395.1.fna.gz\tu.fasta\t96.37\t31.02\t38.48\tShigella sonnei\tcontig_1\n"
                 "GCF_000005845.2.fna.gz\tu.fasta\t95.99\t36.44\t42.02\tEscherichia coli K-12\tcontig_1\n")
    return p


def run_delta(tmp_path, delta, skani):
    out = tmp_path / ("d%s.species.tsv" % delta)
    assert sc.main(["--unit", "u", "--species", "Escherichia coli", "--ani-min", "95", "--skani", str(skani),
                    "--references", str(REFS), "--confusable-delta", str(delta), "--out", str(out)]) == 0
    return list(csv.DictReader(open(out), delimiter="\t"))[0]


def test_close_relative_within_delta_is_unresolved(tmp_path):
    row = run_delta(tmp_path, 0.5, ecoli_shigella(tmp_path))
    assert row["confirmed"] == "false" and row["reason"] == "unresolved_confusable"
    assert row["best_species"] == "Shigella sonnei" and row["panel_ani"] == "95.99"


def test_close_relative_outside_delta_is_other_species(tmp_path):
    row = run_delta(tmp_path, 0.2, ecoli_shigella(tmp_path))
    assert row["reason"] == "best_hit_other_species" and row["panel_ani"] == "95.99"


def test_unresolved_needs_panel_ani_above_min(tmp_path):
    p = tmp_path / "low.skani.tsv"
    p.write_text(SNIPPET.read_text().splitlines()[0] + "\n"
                 "GCF_002950395.1.fna.gz\tu.fasta\t94.30\t31.0\t38.0\tShigella sonnei\tcontig_1\n"
                 "GCF_000005845.2.fna.gz\tu.fasta\t94.00\t36.0\t42.0\tEscherichia coli K-12\tcontig_1\n")
    assert run_delta(tmp_path, 0.5, p)["reason"] == "best_hit_other_species"
