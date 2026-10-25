"""Contract: every amrfinder_organism in the panel is a valid AMRFinderPlus --organism."""
import csv

from conftest import DATA

ASSETS = DATA.parents[1] / "assets"


def test_panel_organisms_are_known_to_amrfinderplus():
    text = (DATA / "pathogen_contracts" / "amrfinder_organisms.txt").read_text()
    known = {o.strip() for o in text.split(":", 1)[-1].replace("\n", ",").split(",") if o.strip()}
    used = {r["amrfinder_organism"] for r in csv.DictReader(open(ASSETS / "pathogen_panel.csv")) if r["amrfinder_organism"]}
    assert used and used <= known, sorted(used - known)
