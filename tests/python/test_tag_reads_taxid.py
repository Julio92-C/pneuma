import pytest
from conftest import DATA, load_bin

tag = load_bin("tag_reads_taxid")


def test_tags_reads_with_taxid():
    with open(DATA / "classified_small.fastq") as fh:
        out = list(tag.tag_records(fh))
    assert out == [
        ">read1_kraken:taxid|1301\nACGTACGT\n",
        ">read2_kraken:taxid|727\nGGCC\n",
    ]


def test_missing_taxid_raises():
    with pytest.raises(ValueError, match="kraken:taxid"):
        list(tag.tag_records(["@r1 nothing\n", "AC\n", "+\n", "II\n"]))


def test_truncated_record_raises():
    with pytest.raises(ValueError, match="truncated"):
        list(tag.tag_records(["@r1 kraken:taxid|1\n", "AC\n"]))


def test_cli_writes_fasta(tmp_path):
    out = tmp_path / "S1_kraken2ID_PFP.fasta"
    tag.main([str(DATA / "classified_small.fastq"), str(out)])
    assert out.read_text().startswith(">read1_kraken:taxid|1301\n")
