import pytest
from conftest import BIN, DATA, load_bin

mm = load_bin("make_meridian_config")
TPL_PATH = BIN.parent / "assets" / "meridian_config.template.yaml"
TEMPLATE = TPL_PATH.read_text(encoding="utf-8")


def test_render_fills_all_placeholders():
    text = mm.render(TEMPLATE, "lung", "metadata_test.csv", "sample_type", "ID", ["CNF0-0001"])
    assert "{{" not in text
    assert 'group_cols:\n    - "sample_type"' in text
    assert 'random_effect: "ID"' in text
    assert 'controls: ["CNF0-0001"]' in text
    assert 'file: "inputs/metadata_test.csv"' in text
    assert 'id: "lung"' in text
    assert 'name: "lung (PNEUMA)"' in text
    assert "__PROJECT_ROOT__" in text


def test_render_defaults():
    text = mm.render(TEMPLATE, "lung", "m.csv", "sample_type", None, [])
    assert "random_effect: ~" in text
    assert "controls: []" in text


def test_render_quotes_yaml_keywords():
    text = mm.render(TEMPLATE, "2024", "m.csv", "yes", None, [])
    assert 'id: "2024"' in text
    assert 'group_cols:\n    - "yes"' in text
    assert 'fixed_effects:\n    - "yes"' in text
    assert "random_effect: ~" in text


def test_render_escapes_quotes_and_backslashes():
    text = mm.render(TEMPLATE, 'a"b: c\\d', "m.csv", "g", None, [])
    assert 'id: "a\\"b: c\\\\d"' in text


def test_rendered_yaml_types():
    yaml = pytest.importorskip("yaml")
    cfg = yaml.safe_load(mm.render(TEMPLATE, "2024", "m.csv", "yes", "no", ["null", "on"]))
    assert cfg["study"]["id"] == "2024"
    assert cfg["metadata"]["group_cols"] == ["yes"]
    assert cfg["metadata"]["fixed_effects"] == ["yes"]
    assert cfg["metadata"]["random_effect"] == "no"
    assert cfg["metadata"]["controls"] == ["null", "on"]
    cfg2 = yaml.safe_load(mm.render(TEMPLATE, 'a"b: c', "m.csv", "g", None, []))
    assert cfg2["study"]["id"] == 'a"b: c'
    assert cfg2["metadata"]["random_effect"] is None


def test_check_metadata_ok():
    mm.check_metadata(DATA / "metadata_test.csv", "sample_type", "ID", ["CNF0-0001"])


def test_check_metadata_bom(tmp_path):
    p = tmp_path / "meta.csv"
    p.write_bytes("﻿sample,sample_type\nA,x\nB,y\nC,z\n".encode("utf-8"))
    mm.check_metadata(p, "sample_type", None, ["A"])


def test_check_metadata_no_sample_column(tmp_path):
    p = tmp_path / "meta.csv"
    p.write_text("id,sample_type\nA,x\n", encoding="utf-8")
    with pytest.raises(ValueError, match="sample"):
        mm.check_metadata(p, "sample_type", None, [])


@pytest.mark.parametrize("group,re_col,controls,msg", [
    ("missing_col", None, [], "group column"),
    ("sample_type", "nope", [], "random-effect column"),
    ("sample_type", None, ["NOT-THERE"], "control"),
])
def test_check_metadata_errors(group, re_col, controls, msg):
    with pytest.raises(ValueError, match=msg):
        mm.check_metadata(DATA / "metadata_test.csv", group, re_col, controls)


def _cli(out, controls):
    mm.main(["--template", str(TPL_PATH),
             "--out", str(out), "--study-id", "lung",
             "--metadata", str(DATA / "metadata_test.csv"),
             "--group", "sample_type", "--controls", controls])


def test_cli(tmp_path):
    out = tmp_path / "config.yaml"
    _cli(out, "CNF0-0001")
    assert 'controls: ["CNF0-0001"]' in out.read_text(encoding="utf-8")


def test_cli_strips_control_whitespace(tmp_path):
    out = tmp_path / "config.yaml"
    _cli(out, "CNF0-0001, ")
    assert 'controls: ["CNF0-0001"]' in out.read_text(encoding="utf-8")


@pytest.mark.parametrize("body,controls", [
    ("sample,g\nA,x\nB,x\nC,x\n", []),
    ("sample,g\nA,x\nB,NA\nC,\n", []),
    ("sample,g\nA,x\nB,x\nCTRL,z\n", ["CTRL"]),   # controls do not count as a level
])
def test_check_metadata_needs_two_levels(tmp_path, body, controls):
    p = tmp_path / "meta.csv"
    p.write_text(body, encoding="utf-8")
    with pytest.raises(ValueError, match="'g'"):
        mm.check_metadata(p, "g", None, controls)


def test_check_metadata_two_levels_ok(tmp_path):
    p = tmp_path / "meta.csv"
    p.write_text("sample,g\nA,x\nB,y\nCTRL,z\n", encoding="utf-8")
    mm.check_metadata(p, "g", None, ["CTRL"])


# MERIDIAN >= 1.2.0 options: control matching, patient-aware tests, beta read floor.
def test_render_v12_defaults_keep_old_behaviour():
    yaml = pytest.importorskip("yaml")
    cfg = yaml.safe_load(mm.render(TEMPLATE, "lung", "m.csv", "sample_type", None, []))
    assert cfg["metadata"]["control_match_cols"] == []
    assert cfg["stats"]["random_effect_tests"] is False
    assert cfg["beta_diversity"]["min_sample_reads"] == 0


def test_render_v12_options():
    yaml = pytest.importorskip("yaml")
    cfg = yaml.safe_load(mm.render(TEMPLATE, "lung", "m.csv", "sample_type", "ID", [],
                                   control_match_cols=["seq_run", "batch"],
                                   random_effect_tests=True, min_sample_reads=100))
    assert cfg["metadata"]["control_match_cols"] == ["seq_run", "batch"]
    assert cfg["stats"]["random_effect_tests"] is True
    assert cfg["beta_diversity"]["min_sample_reads"] == 100


def test_check_metadata_match_col_missing():
    with pytest.raises(ValueError, match="control-match column 'seq_run'"):
        mm.check_metadata(DATA / "metadata_test.csv", "sample_type", "ID", [],
                          control_match_cols=["seq_run"])


def test_cli_v12_options(tmp_path):
    yaml = pytest.importorskip("yaml")
    out = tmp_path / "config.yaml"
    mm.main(["--template", str(TPL_PATH), "--out", str(out), "--study-id", "lung",
             "--metadata", str(DATA / "metadata_test.csv"), "--group", "sample_type",
             "--random-effect", "ID", "--random-effect-tests",
             "--control-match-cols", "ID", "--min-sample-reads", "100"])
    cfg = yaml.safe_load(out.read_text(encoding="utf-8"))
    assert cfg["metadata"]["control_match_cols"] == ["ID"]
    assert cfg["stats"]["random_effect_tests"] is True
    assert cfg["beta_diversity"]["min_sample_reads"] == 100


def test_cli_random_effect_tests_need_random_effect(tmp_path):
    with pytest.raises(ValueError, match="random_effect_tests needs"):
        mm.main(["--template", str(TPL_PATH), "--out", str(tmp_path / "c.yaml"),
                 "--study-id", "lung", "--metadata", str(DATA / "metadata_test.csv"),
                 "--group", "sample_type", "--random-effect-tests"])
