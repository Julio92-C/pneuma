"""resolve_panel.py: names -> taxids (+children), synonyms, suggestions, panel validation."""
import pytest

from conftest import DATA, load_bin

rp = load_bin("resolve_panel")
ASSETS = DATA.parents[1] / "assets"
HEADER = "species,tier,who_category,action,reference_accession,confusable_with,ani_min,amrfinder_organism\n"


def taxdump(tmp_path):
    d = tmp_path / "taxdump"
    d.mkdir()
    nodes = [("1", "1", "no rank"), ("1301", "1", "genus"), ("1313", "1301", "species"),
             ("170187", "1313", "no rank"), ("28037", "1301", "species"),
             ("3705180", "1", "genus"), ("729", "3705180", "species")]
    names = [("1", "root", "scientific name"), ("1301", "Streptococcus", "scientific name"),
             ("1313", "Streptococcus pneumoniae", "scientific name"),
             ("170187", "Streptococcus pneumoniae TIGR4", "scientific name"),
             ("28037", "Streptococcus mitis", "scientific name"),
             ("3705180", "Parahaemophilus", "scientific name"),
             ("729", "Parahaemophilus parainfluenzae", "scientific name"),
             ("729", "Haemophilus parainfluenzae", "synonym")]
    (d / "nodes.dmp").write_text("".join("%s\t|\t%s\t|\t%s\t|\t\t|\n" % n for n in nodes))
    (d / "names.dmp").write_text("".join("%s\t|\t%s\t|\t\t|\t%s\t|\n" % n for n in names))
    return d


def refs(tmp_path, rows=(("Streptococcus pneumoniae", "GCF_1"), ("Streptococcus mitis", "GCF_2"),
                         ("Haemophilus parainfluenzae", "GCF_3"))):
    p = tmp_path / "refs.tsv"
    p.write_text("species\taccession\trole\n" + "".join("%s\t%s\tpanel\n" % r for r in rows))
    return p


def panel(tmp_path, body):
    p = tmp_path / "panel.csv"
    p.write_text(HEADER + body)
    return p


def test_resolves_children_synonym_and_confusables(tmp_path):
    out = tmp_path / "out.tsv"
    rows = rp.run(panel(tmp_path,
                        "Streptococcus pneumoniae,1,BPPL 2024 medium,reconstruct,GCF_1,Streptococcus mitis,,Streptococcus_pneumoniae\n"
                        "Haemophilus parainfluenzae,R,,reconstruct,GCF_3,,97,\n"),
                  refs(tmp_path), taxdump(tmp_path), out)
    spn, hpa = rows
    assert spn["taxid"] == "1313" and spn["taxids"] == "1313;170187"
    assert spn["slug"] == "Streptococcus_pneumoniae" and spn["ani_min"] == "95.0"
    assert spn["confusable_accessions"] == "GCF_2"
    assert hpa["taxid"] == "729" and hpa["matched_name_class"] == "synonym"
    assert hpa["scientific_name"] == "Parahaemophilus parainfluenzae" and hpa["ani_min"] == "97.0"
    assert out.read_text().splitlines()[0].split("\t") == rp.OUT_COLUMNS


def test_misspelled_name_fails_with_suggestion(tmp_path):
    with pytest.raises(ValueError) as exc:
        rp.run(panel(tmp_path, "Streptococus pneumoniae,1,,reconstruct,GCF_1,,95,\n"),
               refs(tmp_path, (("Streptococus pneumoniae", "GCF_1"),)), taxdump(tmp_path), tmp_path / "o.tsv")
    assert "not found" in str(exc.value) and "did you mean: Streptococcus pneumoniae" in str(exc.value)


def test_non_species_rank_and_bad_columns_fail_together(tmp_path):
    with pytest.raises(ValueError) as exc:
        rp.run(panel(tmp_path,
                     "Streptococcus,1,,reconstruct,GCF_9,Streptococcus nonexistentii,95,\n"
                     "Streptococcus mitis,X,,assemble,,,abc,\n"),
               refs(tmp_path), taxdump(tmp_path), tmp_path / "o.tsv")
    msg = str(exc.value)
    assert "rank 'genus', expected species" in msg
    assert "does not match" in msg                   # GCF_9 vs references (absent)
    assert "confusable species without an accession" in msg
    assert "tier 'X'" in msg and "action 'assemble'" in msg


def test_reconstruct_needs_reference(tmp_path):
    with pytest.raises(ValueError, match="reconstruct needs reference_accession"):
        rp.run(panel(tmp_path, "Streptococcus mitis,R,,reconstruct,,,95,\n"),
               refs(tmp_path), taxdump(tmp_path), tmp_path / "o.tsv")


def test_asset_panel_resolves_against_fixture_taxdump(tmp_path):
    rows = rp.run(ASSETS / "pathogen_panel.csv", ASSETS / "panel_references.tsv",
                  DATA / "pathogen_fixture" / "taxdump", tmp_path / "o.tsv")
    assert len(rows) == 25
    assert sum(r["action"] == "detect" for r in rows) == 5
    assert all(r["tier"] in rp.TIERS for r in rows)
    assert {r["species"]: r["taxid"] for r in rows}["Haemophilus parainfluenzae"] == "729"


def test_cli_exit_code(tmp_path):
    rc = rp.main(["--panel", str(DATA / "pathogen_contracts" / "panel_misspelled.csv"),
                  "--references", str(ASSETS / "panel_references.tsv"),
                  "--taxdump", str(DATA / "pathogen_fixture" / "taxdump"), "--out", str(tmp_path / "o.tsv")])
    assert rc == 1
