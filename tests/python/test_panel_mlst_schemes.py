"""Contract: every mlst_scheme in the panel is a scheme of the pinned mlst container (or 'none')."""
import csv

from conftest import DATA

ASSETS = DATA.parents[1] / "assets"


def test_panel_mlst_schemes_exist_in_container():
    # captured with: mlst --list  (quay.io/biocontainers/mlst:2.25.0--hdfd78af_0)
    known = set((DATA / "pathogen_contracts" / "mlst_schemes.txt").read_text().split())
    used = {r["mlst_scheme"] for r in csv.DictReader(open(ASSETS / "pathogen_panel.csv")) if r["mlst_scheme"]}
    assert used, "panel has no mlst_scheme values"
    assert used - {"none"} <= known, sorted(used - {"none"} - known)
