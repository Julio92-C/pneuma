"""The public-data pathogen fixture must stay reproducible, realistic and complete."""
import csv
import gzip
import importlib.util
import re
from pathlib import Path

from conftest import DATA

FIX = DATA / "pathogen_fixture"
MINI = DATA / "pathogen_mini"
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9-]*$")


def _gen():
    spec = importlib.util.spec_from_file_location(
        "make_pathogen_fixture", Path(__file__).resolve().parents[1] / "make_pathogen_fixture.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_taxdump_is_rebuilt_identically_from_lineages():
    gen = _gen()
    nodes, names = gen.build_taxdump(gen.load_lineages(), gen.BASE_TAXDUMP)
    assert "\n".join(nodes) + "\n" == (FIX / "taxdump" / "nodes.dmp").read_text()
    assert "\n".join(names) + "\n" == (FIX / "taxdump" / "names.dmp").read_text()


def test_taxdump_keeps_the_base_and_has_the_renamed_species():
    names = (FIX / "taxdump" / "names.dmp").read_text()
    for base_line in (DATA / "taxdump" / "names.dmp").read_text().splitlines():
        assert base_line in names            # superset of the Phase 1 test taxdump (Recentrifuge)
    assert "729\t|\tHaemophilus parainfluenzae\t|\t\t|\tsynonym\t|" in names
    assert "\n170187\t|\t1313\t|\tno rank\t|" in "\n" + (FIX / "taxdump" / "nodes.dmp").read_text()


def test_kreport_roundtrip():
    gen = _gen()
    parent, rank, name = gen.read_taxdump(FIX / "taxdump")
    lines = gen.kreport({"1313": 7, "170187": 3, "1301": 2, "U": 5}, parent, rank, name)
    rows = [line.split("\t") for line in lines]
    assert rows[0][1] == "5" and rows[0][3:6] == ["U", "0", "unclassified"]
    by_taxid = {r[4]: r for r in rows}
    assert by_taxid["1"][1] == "12"                    # root clade = classified reads
    assert by_taxid["1313"][1] == "10" and by_taxid["1313"][2] == "7"
    assert by_taxid["170187"][3] == "S1"               # strain below species
    assert by_taxid["131567"][3] == "R1"
    assert by_taxid["1313"][5].strip() == "Streptococcus pneumoniae"
    assert by_taxid["1313"][5].startswith("  " * 9)    # root, cellular org., D, K, P, C, O, F, G


def test_design_covers_every_required_scenario():
    rows = list(csv.DictReader(open(FIX / "expected_status.tsv"), delimiter="\t"))
    assert {"TYPED", "LOW_COVERAGE", "SPECIES_NOT_CONFIRMED", "MIXED_STRAINS"} <= {r["status"] for r in rows}
    patient = [r for r in rows if r["level"] == "patient"]
    assert any(r["status"] == "TYPED" and r["contributing_samples"] == "PD1-S1;PD1-S2" for r in patient)
    assert any(r["status"] == "MIXED_STRAINS" and r["contributing_samples"] == "PE1-S1;PE1-S2" for r in patient)
    gen = _gen()
    assert any(control for _s, _p, control, _c in gen.DESIGN)   # control: NOT_DETECTED everywhere
    for sample, patient_id, _control, _components in gen.DESIGN:
        assert ID_RE.match(sample)
        assert patient_id == "" or ID_RE.match(patient_id)


def test_mini_fixture_is_consistent():
    gen = _gen()
    samples = list(csv.DictReader(open(MINI / "samples.tsv"), delimiter="\t"))
    assert [s["sample"] for s in samples] == ["MINI-S1", "MINI-S2", "MINI-S3"]
    for s in samples:
        sid = s["sample"]
        with gzip.open(MINI / f"{sid}.classified.fastq.gz", "rt") as fh:
            fq_ids = [line[1:].split()[0] for i, line in enumerate(fh) if i % 4 == 0]
        assignments = (MINI / f"{sid}.kraken2.classifiedreads.txt").read_text().splitlines()
        assert fq_ids == [a.split("\t")[1] for a in assignments if a.startswith("C\t")]
        assert int(s["n_reads"]) == len(assignments)
        assert (MINI / "bracken" / f"{sid}.tsv").read_text().startswith(gen.BRACKEN_HEADER)
    assert (MINI / "ref" / "MINI-TIGR4-SEG.fna.gz").exists()
    index = (MINI / "ref" / "panel_references_index.tsv").read_text().splitlines()
    assert index[1].split("\t")[1] == "MINI-TIGR4-SEG"
