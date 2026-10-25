"""pathogen_report.py: exactly one status per sample x species, patient-level rows, HTML."""
import csv
import filecmp
import importlib.util
from pathlib import Path

import pytest

from conftest import DATA, load_bin

pr = load_bin("pathogen_report")
spec = importlib.util.spec_from_file_location("report_fixture", Path(__file__).with_name("report_fixture.py"))
rf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rf)
SPN, HIN, MTB, SMI = rf.SPN, rf.HIN, rf.MTB, rf.SMI


def build(tmp_path):
    root = rf.write(tmp_path / "in")
    out = tmp_path / "out"
    rc = pr.main(["--panel", str(root / "panel_resolved.tsv"), "--candidates", str(root / "pathogen_candidates.tsv"),
                  "--samples", str(root / "pathogen_samples.tsv"), "--units", str(root / "units.tsv"),
                  "--results", str(root), "--panel-file", str(root / "pathogen_panel.csv"),
                  "--manifest", str(root / "database_manifest.tsv"), "--pipeline-version", "0.3.0",
                  "--gate-min-depth", "10", "--gate-min-breadth", "0.5", "--outdir", str(out)])
    assert rc == 0
    rows = list(csv.DictReader(open(out / "pathogen_cohort.tsv"), delimiter="\t"))
    return out, rows


def status(rows, sample, species):
    hits = [r for r in rows if r["level"] == "sample" and r["sample"] == sample and r["species"] == species]
    assert len(hits) == 1
    return hits[0]


def test_one_status_per_sample_and_species(tmp_path):
    _out, rows = build(tmp_path)
    sample_rows = [r for r in rows if r["level"] == "sample"]
    assert len(sample_rows) == len(rf.SAMPLES) * len(rf.PANEL)
    assert len({(r["sample"], r["species"]) for r in sample_rows}) == len(sample_rows)
    assert all(r["status"] in pr.STATUSES for r in rows)
    assert list(rows[0]) == pr.COLUMNS


def test_every_status(tmp_path):
    _out, rows = build(tmp_path)
    s1 = status(rows, "S1", SPN)
    assert (s1["status"], s1["st"], s1["scheme"]) == ("TYPED", "1", "spneumoniae")
    assert s1["amr_genes"] == "tet(M)" and s1["vf_genes"] == "ply" and s1["gene_source"] == "contig"
    assert s1["completeness"] == "97.5" and s1["n50"] == "900000"
    assert status(rows, "S1", HIN)["status"] == "NOT_DETECTED"
    assert status(rows, "S1", MTB)["status"] == "DETECTED_ONLY"
    s2 = status(rows, "S2", SPN)
    assert s2["status"] == "LOW_COVERAGE" and s2["amr_genes"] == "tet(M)" and s2["gene_source"] == "read-level"
    assert "not extracted" in s2["notes"]
    assert status(rows, "S2", SMI)["status"] == "ASSEMBLY_FAILED"
    s3 = status(rows, "S3", SPN)
    assert s3["status"] == "SPECIES_NOT_CONFIRMED" and s3["st"] == "" and s3["amr_genes"] == ""
    assert s3["ani_reference"] == "Streptococcus mitis (GCF_018603625.2)"
    assert status(rows, "S3", HIN)["status"] == "PARTIAL"
    assert status(rows, "S3", SMI)["status"] == "NO_SCHEME"
    assert all(status(rows, "S8", sp)["status"] == "LOW_DEPTH" for sp in (SPN, HIN, MTB, SMI))
    assert all(status(rows, "CTRL", sp)["status"] == "NOT_DETECTED" for sp in (SPN, HIN, MTB, SMI))


def test_patient_level_rows_and_pooled_references(tmp_path):
    _out, rows = build(tmp_path)
    patient = {(r["patient"], r["species"]): r for r in rows if r["level"] == "patient"}
    p4 = patient[("P4", HIN)]
    assert p4["status"] == "TYPED" and p4["st"] == "3" and p4["contributing_samples"] == "S4;S5"
    assert p4["sample"] == "" and p4["reads"] == "780"
    p5 = patient[("P5", SPN)]
    assert p5["status"] == "MIXED_STRAINS" and p5["st"] == "" and p5["mixed_het_per_kb"] == "5.000"
    for s in ("S4", "S5"):
        row = status(rows, s, HIN)
        assert row["status"] == "LOW_COVERAGE" and row["pooled_result"] == "TYPED (pooled, patient P4)"
        assert row["st"] == ""                       # never presented as this sample's genome
    assert status(rows, "S6", SPN)["pooled_result"] == "MIXED_STRAINS (pooled, patient P5)"


def test_html_per_sample_with_patient_section(tmp_path):
    out, _rows = build(tmp_path)
    pages = sorted(p.name for p in (out / "html").iterdir())
    assert pages == sorted("%s.pathogens.html" % s[0] for s in rf.SAMPLES)
    s4 = (out / "html" / "S4.pathogens.html").read_text()
    assert "Patient-level (pooled) results, patient P4" in s4 and "not genomes from this sample alone" in s4
    for section in ("Tier 1-3 priority pathogens", "Detect-only (Tier D)", "Tier R commensal AMR/VF reservoirs"):
        assert section in s4
    assert "Patient-level" not in (out / "html" / "S1.pathogens.html").read_text()


def test_missing_gate_is_a_wiring_error(tmp_path):
    root = rf.write(tmp_path / "in")
    (root / "gate" / "S1.Streptococcus_pneumoniae.gate.tsv").unlink()
    with pytest.raises(pr.ContractError, match="no gate result for unit S1.Streptococcus_pneumoniae"):
        pr.build_report(root / "panel_resolved.tsv", root / "pathogen_candidates.tsv", root / "pathogen_samples.tsv",
                        root / "units.tsv", root)


def test_committed_nf_test_fixture_is_current(tmp_path):
    rf.write(tmp_path / "fx")
    committed = DATA / "pathogen_report"
    for p in sorted((tmp_path / "fx").rglob("*")):
        if p.is_file():
            rel = p.relative_to(tmp_path / "fx")
            assert filecmp.cmp(p, committed / rel, shallow=False), "regenerate: python3 tests/python/report_fixture.py tests/data/pathogen_report (%s)" % rel


# ---- Task 10 additions: real-contract parsers, missing optional outputs ----

def test_parsers_match_real_captured_headers():
    c = DATA / "pathogen_contracts"
    assert pr.parse_abricate_tab(c / "abricate_contigs.txt") == ["zmpC"]   # data row has 14 of 15 fields
    assert pr.parse_amrfinder(c / "amrfinderplus.tsv") == ["tet(M)"]
    assert pr.parse_mlst(c / "mlst.tsv")[:3] == ("TYPED", "spneumoniae", "205")
    assert pr.parse_checkm2([c / "checkm2_report.tsv"])["spn_tigr4"] == ("99.9", "0.3")
    assert pr.parse_quast([c / "quast_report.tsv"])["MINI-S1.Streptococcus_pneumoniae"] == {
        "contigs": "1", "total_length": "122328", "n50": "122328"}


def test_missing_mobsuite_output_means_no_plasmids(tmp_path):
    root = rf.write(tmp_path / "in")
    assert not (root / "mobsuite").exists()          # MOB-suite writes nothing for plasmid-free assemblies
    assert status(build_rows(root), "S1", SPN)["mge_genes"] == ""


def test_mobsuite_replicons_are_reported_and_short_rows_tolerated(tmp_path):
    root = rf.write(tmp_path / "in")
    (root / "mobsuite").mkdir()
    (root / "mobsuite" / "S1.Streptococcus_pneumoniae.mobtyper_results.txt").write_text(
        "sample_id\tnum_contigs\tsize\trep_type(s)\trep_type_accession(s)\n"
        "contig_1\t1\t5000\tIncFIB,rep1\n")                # trailing column stripped
    assert status(build_rows(root), "S1", SPN)["mge_genes"] == "mob:IncFIB;mob:rep1"


def build_rows(root):
    rows, _p, _s = pr.build_report(root / "panel_resolved.tsv", root / "pathogen_candidates.tsv",
                                   root / "pathogen_samples.tsv", root / "units.tsv", root)
    return rows
