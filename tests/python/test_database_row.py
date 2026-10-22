"""database_row.py: one manifest row with a stable sha256 for files and directories."""
import json
import tarfile

import pytest

from conftest import load_bin

dr = load_bin("database_row")


def test_directory_digest_is_stable_and_content_based(tmp_path):
    a = tmp_path / "a"
    (a / "sub").mkdir(parents=True)
    (a / "x.txt").write_text("one")
    (a / "sub" / "y.txt").write_text("two")
    first = dr.sha256_path(a)
    (a / "x.txt").touch()                         # timestamp change only
    assert dr.sha256_path(a) == first
    (a / "x.txt").write_text("changed")
    assert dr.sha256_path(a) != first


def test_bakta_and_amrfinderplus_versions(tmp_path):
    db = tmp_path / "db-light"
    db.mkdir()
    (db / "version.json").write_text(json.dumps({"date": "2025-02-24", "major": 6, "minor": 0, "type": "light"}))
    assert dr.version_from("bakta", db) == "6.0-light (2025-02-24)"
    src = tmp_path / "version.txt"
    src.write_text("2026-03-24.1\n")
    tgz = tmp_path / "amrfinderdb.tar.gz"
    with tarfile.open(tgz, "w:gz") as t:
        t.add(src, arcname="./version.txt")
    assert dr.version_from("amrfinderplus", tgz) == "2026-03-24.1"


def test_cli_writes_one_row_and_rejects_empty_version(tmp_path):
    f = tmp_path / "m.dmnd"
    f.write_text("x")
    out = tmp_path / "07_checkm2_db.row.tsv"
    assert dr.main(["--name", "checkm2_db", "--version", "3", "--rel", "checkm2/m.dmnd",
                    "--target", str(f), "--out", str(out)]) == 0
    name, version, rel, sha = out.read_text().rstrip("\n").split("\t")
    assert (name, version, rel) == ("checkm2_db", "3", "checkm2/m.dmnd") and len(sha) == 64
    with pytest.raises(SystemExit):
        dr.main(["--name", "x", "--version", "", "--rel", "r", "--target", str(f), "--out", str(out)])
