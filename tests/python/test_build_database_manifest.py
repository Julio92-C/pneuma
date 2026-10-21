import subprocess
import sys

import pytest
from conftest import BIN, load_bin

mod = load_bin("build_database_manifest")
SHA = "a" * 64


def write(tmp_path, name, *rows):
    f = tmp_path / name
    f.write_text("".join("\t".join(r) + "\n" for r in rows))
    return str(f)


def run_cli(root, out, *files):
    return subprocess.run([sys.executable, str(BIN / "build_database_manifest.py"), "--root", root,
                           "--out", str(out), *files], capture_output=True, text=True)


def test_orders_rows_and_prefixes_root(tmp_path):
    ab = write(tmp_path, "04.tsv", ("abricate_vfdb", "2025-Jan-01", "abricate/vfdb", SHA),
               ("abricate_card", "2025-Jan-01", "abricate/card", SHA))
    ho = write(tmp_path, "03.tsv", ("hostile_index", "human-t2t-hla", "hostile/x.mmi", SHA))
    tx = write(tmp_path, "02.tsv", ("ncbi_taxdump", "2025-04", "taxdump_2025-04", SHA))
    kr = write(tmp_path, "01.tsv", ("kraken2_pluspf", "20250402", "kraken2_pluspf_20250402", SHA))
    rows = mod.build([ab, ho, tx, kr], "/data/db/")
    assert [r[0] for r in rows] == ["kraken2_pluspf", "ncbi_taxdump", "hostile_index",
                                    "abricate_card", "abricate_vfdb"]
    assert rows[0][2] == "/data/db/kraken2_pluspf_20250402"


def test_cli_writes_header_and_tsv(tmp_path):
    kr = write(tmp_path, "01.tsv", ("kraken2_pluspf", "20250402", "k", SHA))
    out = tmp_path / "database_manifest.tsv"
    r = run_cli("/r", out, kr)
    assert r.returncode == 0, r.stderr
    assert out.read_text() == f"name\tversion\tpath\tsha256\nkraken2_pluspf\t20250402\t/r/k\t{SHA}\n"


@pytest.mark.parametrize("row", [
    ("kraken2_pluspf", "20250402", "k", ""),            # empty sha
    ("abricate_card", "", "abricate/card", SHA),        # empty date
    ("abricate_card", "2025-Jan-01", "abricate/card"),  # missing column
])
def test_incomplete_row_fails(tmp_path, row):
    f = write(tmp_path, "04.tsv", row)
    out = tmp_path / "m.tsv"
    r = run_cli("/r", out, f)
    assert r.returncode == 1
    assert "incomplete row" in r.stderr
    assert not out.exists()


def test_empty_input_fails(tmp_path):
    f = tmp_path / "empty.tsv"
    f.write_text("")
    with pytest.raises(ValueError, match="no database rows"):
        mod.build([str(f)], "/r")
