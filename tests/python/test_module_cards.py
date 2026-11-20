"""docs/images/make_module_cards.py: one SVG card per subworkflow, valid, readable and up to date."""
import importlib.util
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("make_module_cards", ROOT / "docs" / "images" / "make_module_cards.py")
mc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mc)
SVG = "{http://www.w3.org/2000/svg}"
SLUGS = ["preprocess", "taxonomy", "contamination", "functional", "genomic_context", "mobilome", "pathogens", "meridian"]


def test_one_card_per_subworkflow():
    assert [c["slug"] for c in mc.CARDS] == SLUGS


def test_cards_are_valid_svg_with_every_step_input_and_output():
    for card in mc.CARDS:
        root = ET.fromstring(mc.render(card))
        texts = " ".join(t.text or "" for t in root.iter(SVG + "text"))
        assert card["title"] in texts
        for s in card["steps"]:
            assert s["tool"] in texts and s["question"] in texts
        for _kind, label in card["inputs"] + card["outputs"]:
            assert label in texts


def test_optional_steps_are_dashed():
    for card in mc.CARDS:
        assert mc.render(card).count('stroke-dasharray="7 5"') == sum(s["optional"] for s in card["steps"])
    assert any(s["optional"] for c in mc.CARDS for s in c["steps"])


def test_labels_fit_their_boxes():
    for card in mc.CARDS:
        assert all(kind in ("reads", "seqs", "table", "db") for kind, _l in card["inputs"] + card["outputs"])
        assert all(len(label) <= 32 for _k, label in card["inputs"] + card["outputs"])
        assert all(len(s["tool"]) <= 40 and len(s["question"]) <= 70 and len(s["blind"]) <= 76 for s in card["steps"])


def test_committed_cards_are_current(tmp_path):
    assert mc.main(["--outdir", str(tmp_path)]) == 0
    for slug in SLUGS:
        new = (tmp_path / (slug + ".svg")).read_text(encoding="utf-8")
        assert (ROOT / "docs" / "images" / "cards" / (slug + ".svg")).read_text(encoding="utf-8") == new
