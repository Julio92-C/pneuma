import json
import re
from conftest import BIN

ROOT = BIN.parent


def test_abricate_db_param_declared():
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    props = {}
    for grp in schema.get("$defs", schema.get("definitions", {})).values():
        props.update(grp.get("properties", {}))
    props.update(schema.get("properties", {}))
    assert "abricate_db" in props
    assert props["abricate_db"]["type"] == "string"
    cfg = (ROOT / "nextflow.config").read_text()
    assert re.search(r"^\s*abricate_db\s*=\s*null\b", cfg, re.M)
