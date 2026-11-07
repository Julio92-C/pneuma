"""hit_reads.py: per-sample hit reads, hits table, lengths and taxa for the genomic-context analysis."""
from conftest import load_bin

hr = load_bin("hit_reads")
HEAD = "#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE\n"


def test_hit_reads(tmp_path):
    fa = tmp_path / "S1_kraken2ID_PFP.fasta"
    fa.write_text(">r1_kraken:taxid|562\nACGTACGTAC\n>r2_kraken:taxid|1280\nGGGG\n>r3_kraken:taxid|562\nTTTTTT\n")
    csvf = tmp_path / "S1_abricate.csv"
    csvf.write_text(HEAD + "S1_kraken2ID_PFP.fasta,r1_kraken:taxid|562,2,8,+,tet(M),1-7/7,===,0/0,100.00,99.00,card,X,p,tetracycline\n"
                    + HEAD + "S1_kraken2ID_PFP.fasta,r3_kraken:taxid|562,1,6,+,ColRNAI_1,1-6/6,===,0/0,100.00,99.00,plasmidfinder,Y,p,\n")
    assert hr.main(["--sample", "S1", "--fasta", str(fa), "--abricate", str(csvf), "--out-prefix", str(tmp_path / "S1")]) == 0
    assert (tmp_path / "S1.hit_reads.fa").read_text() == ">S1__r1\nACGTACGTAC\n>S1__r3\nTTTTTT\n"
    assert (tmp_path / "S1.lengths.tsv").read_text() == "S1__r1\t10\nS1__r3\t6\n"
    assert (tmp_path / "S1.taxa.tsv").read_text() == "S1__r1\t562\nS1__r3\t562\n"
    lines = (tmp_path / "S1.hits.tsv").read_text().splitlines()
    assert lines[0].startswith("#FILE\tSEQUENCE\tSTART") and len(lines) == 3
    assert lines[1].split("\t")[1] == "S1__r1" and lines[2].split("\t")[11] == "plasmidfinder"


def test_no_hits(tmp_path):
    fa = tmp_path / "S2.fasta"
    fa.write_text(">r1_kraken:taxid|562\nACGT\n")
    csvf = tmp_path / "S2.csv"
    csvf.write_text(HEAD)
    assert hr.main(["--sample", "S2", "--fasta", str(fa), "--abricate", str(csvf), "--out-prefix", str(tmp_path / "S2")]) == 0
    assert (tmp_path / "S2.hit_reads.fa").read_text() == ""
    assert (tmp_path / "S2.taxa.tsv").read_text() == ""
    assert (tmp_path / "S2.hits.tsv").read_text().count("\n") == 1
