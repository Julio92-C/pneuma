import pytest
from conftest import BIN, DATA, load_bin

mm = load_bin("make_meridian_config")
TEMPLATE = (BIN.parent / "assets" / "meridian_config.template.yaml").read_text()


def test_render_fills_all_placeholders():
    text = mm.render(TEMPLATE, "lung", "metadata_test.csv", "sample_type", "ID", ["CNF0-0001"])
    assert "{{" not in text
    assert "group_cols:\n    - sample_type" in text
    assert "random_effect: ID" in text
    assert 'controls: ["CNF0-0001"]' in text
    assert 'file: "inputs/metadata_test.csv"' in text
    assert "__PROJECT_ROOT__" in text


def test_render_defaults():
    text = mm.render(TEMPLATE, "lung", "m.csv", "sample_type", None, [])
    assert "random_effect: ~" in text
    assert "controls: []" in text


def test_check_metadata_ok():
    mm.check_metadata(DATA / "metadata_test.csv", "sample_type", "ID", ["CNF0-0001"])


@pytest.mark.parametrize("group,re_col,controls,msg", [
    ("missing_col", None, [], "group column"),
    ("sample_type", "nope", [], "random-effect column"),
    ("sample_type", None, ["NOT-THERE"], "control"),
])
def test_check_metadata_errors(group, re_col, controls, msg):
    with pytest.raises(ValueError, match=msg):
        mm.check_metadata(DATA / "metadata_test.csv", group, re_col, controls)


def test_cli(tmp_path):
    out = tmp_path / "config.yaml"
    mm.main(["--template", str(BIN.parent / "assets" / "meridian_config.template.yaml"),
             "--out", str(out), "--study-id", "lung",
             "--metadata", str(DATA / "metadata_test.csv"),
             "--group", "sample_type", "--controls", "CNF0-0001"])
    assert 'controls: ["CNF0-0001"]' in out.read_text()
