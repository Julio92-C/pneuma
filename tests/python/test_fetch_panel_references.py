"""fetch_panel_references.py: NCBI Datasets zip -> <accession>.fna.gz + index (network mocked)."""
import gzip
import io
import zipfile

import pytest

from conftest import load_bin

fp = load_bin("fetch_panel_references")


def fake_zip(accession, seq="ACGT"):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("ncbi_dataset/data/%s/%s_ASM_genomic.fna" % (accession, accession), ">c1\n%s\n" % seq)
        z.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}")
    return buf.getvalue()


def test_writes_genomes_and_index(tmp_path):
    refs = tmp_path / "panel_references.tsv"
    refs.write_text("species\taccession\trole\nA a\tGCF_1.1\tpanel\nB b\tGCF_2.1\tconfusable\n")
    out = tmp_path / "panel_references"
    rows = fp.run(refs, out, getter=lambda acc: fake_zip(acc))
    assert [r["accession"] for r in rows] == ["GCF_1.1", "GCF_2.1"]
    with gzip.open(out / "GCF_1.1.fna.gz", "rt") as fh:
        assert fh.read() == ">c1\nACGT\n"
    index = (out / "panel_references_index.tsv").read_text().splitlines()
    assert index[0] == "species\taccession\trole\tfile\tsha256" and len(index) == 3


def test_package_without_fasta_fails(tmp_path):
    refs = tmp_path / "r.tsv"
    refs.write_text("species\taccession\trole\nA a\tGCF_9.1\tpanel\n")
    empty = io.BytesIO()
    with zipfile.ZipFile(empty, "w") as z:
        z.writestr("README.md", "x")
    with pytest.raises(ValueError, match="GCF_9.1"):
        fp.run(refs, tmp_path / "o", getter=lambda acc: empty.getvalue())
