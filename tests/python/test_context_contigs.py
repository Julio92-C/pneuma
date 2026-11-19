"""context_contigs.py: one sample's Flye contigs renamed <sample>__contig_N and tabulated."""
import csv
import gzip

from conftest import load_bin

cc = load_bin("context_contigs")

INFO = ("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
        "contig_1\t12\t30\tY\tN\t1\t*\t1\n"
        "contig_2\t8\t5\tN\tN\t1\t*\t2\n")


def test_rename_and_table(tmp_path):
    fa = tmp_path / "S1.assembly.fasta.gz"
    with gzip.open(fa, "wt") as fh:
        fh.write(">contig_1\nACGTACGT\nACGT\n>contig_2\nGGGGCCCC\n")
    info = tmp_path / "S1.assembly_info.txt"
    info.write_text(INFO)
    out = tmp_path / "S1"
    assert cc.main(["--sample", "S1", "--fasta", str(fa), "--info", str(info), "--out-prefix", str(out)]) == 0
    assert (tmp_path / "S1.contigs.fa").read_text() == ">S1__contig_1\nACGTACGTACGT\n>S1__contig_2\nGGGGCCCC\n"
    rows = list(csv.DictReader(open(tmp_path / "S1.assembly_info.tsv"), delimiter="\t"))
    assert [(r["contig"], r["length"], r["coverage"], r["circular"]) for r in rows] == [
        ("S1__contig_1", "12", "30", "yes"), ("S1__contig_2", "8", "5", "no")]
    assert {r["sample"] for r in rows} == {"S1"}


def test_empty_assembly_gives_header_only(tmp_path):
    fa, info = tmp_path / "e.fasta", tmp_path / "e.txt"
    fa.write_text("")
    info.write_text("")
    assert cc.main(["--sample", "S1", "--fasta", str(fa), "--info", str(info), "--out-prefix", str(tmp_path / "S1")]) == 0
    assert (tmp_path / "S1.contigs.fa").read_text() == ""
    assert (tmp_path / "S1.assembly_info.tsv").read_text() == "\t".join(cc.INFO_COLUMNS) + "\n"
