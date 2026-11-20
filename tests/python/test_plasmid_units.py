"""plasmid_units.py: plasmid drafts grouped across samples by MinHash identity."""
import csv
import random

from conftest import TEST_TAXA, load_bin, write_taxdump

pu = load_bin("plasmid_units")

CAND_HEAD = ("contig\tsample\tlength\tcoverage\trepeat\treads\tplasmid_score\tchromosome_score\thallmarks\tconjugation"
             "\treplicons\tmob_rep\trelaxase\tmpf\tmobility\tARGs\tVFs\tread_taxa\n")


def seq(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_identity_is_rotation_and_strand_invariant():
    a = seq(20000, 1)
    rotated = a[7000:] + a[:7000]
    assert pu.identity(pu.sketch(a), pu.sketch(revcomp(rotated))) > 0.99
    assert pu.identity(pu.sketch(a), pu.sketch(seq(20000, 2))) < 0.5


def test_clusters_respect_identity_and_length_ratio():
    a, b = seq(20000, 1), seq(20000, 2)
    groups = pu.clusters([("S1__c1", a), ("S2__c4", a[5000:] + a[:5000]), ("S1__c2", b), ("S3__c1", a[:12000])])
    assert sorted(sorted(g) for g in groups) == [["S1__c1", "S2__c4"], ["S1__c2"], ["S3__c1"]]


def test_cli_units_and_members(tmp_path):
    a, b = seq(20000, 1), seq(15000, 2)
    fasta = tmp_path / "c.fa"
    fasta.write_text(">S1__contig_1 x\n%s\n>S2__contig_3\n%s\n>S2__contig_5\n%s\n" % (a, a[3000:] + a[:3000], b))
    cand = tmp_path / "c.tsv"
    cand.write_text(CAND_HEAD
                    + "S1__contig_1\tS1\t20000\t40\tno\t10\t0.9\t0.1\t1\t\trep27\t\tMOBV\t\tmobilizable\ttet(L)\t\t1350:6\n"
                    + "S2__contig_3\tS2\t20000\t20\tyes\t5\t0.9\t0.1\t1\t\trep27\t\tMOBV\t\tmobilizable\ttet(L)\t\t1350:2,1578:1\n"
                    + "S2__contig_5\tS2\t15000\t9\tno\t2\t0.8\t0.1\t0\t\t\t\t\t\tnon-mobilizable\tErmB\t\t\n")
    rec = tmp_path / "S1.r1.recruit.tsv"
    rec.write_text("sample\tround\trecruited\ttotal_reads\nS1\t1\t50\t1000000\n")
    rec2 = tmp_path / "S2.r1.recruit.tsv"
    rec2.write_text("sample\tround\trecruited\ttotal_reads\nS2\t1\t20\t500000\n")
    ise = tmp_path / "i.tsv"
    ise.write_text("seqID\tfamily\tisBegin\tisEnd\tstrand\nS1__contig_1\tIS6\t10\t900\t+\n")
    prefix = tmp_path / "plasmid_units"
    assert pu.main(["--candidates", str(cand), "--fasta", str(fasta), "--isescan", str(ise),
                    "--recruit", str(rec), str(rec2), "--out-prefix", str(prefix)]) == 0
    units = list(csv.DictReader(open(str(prefix) + ".tsv"), delimiter="\t"))
    assert [u["unit"] for u in units] == ["PU001", "PU002"]
    u = units[0]
    assert (u["n_samples"], u["samples"], u["n_contigs"], u["representative"]) == ("2", "S1,S2", "2", "S1__contig_1")
    assert (u["replicons"], u["relaxase"], u["mobility"], u["ARGs"], u["IS_families"]) == (
        "rep27", "MOBV", "mobilizable", "tet(L)", "IS6")
    assert u["read_taxa"] == "1350:8,1578:1"
    assert u["mean_rpm"] == "10.0"                                   # S1 10/1e6*1e6 = 10, S2 5/5e5*1e6 = 10
    members = list(csv.DictReader(open(str(prefix) + "_members.tsv"), delimiter="\t"))
    assert {(m["contig"], m["unit"], m["rpm"]) for m in members} == {
        ("S1__contig_1", "PU001", "10.0"), ("S2__contig_3", "PU001", "10.0"), ("S2__contig_5", "PU002", "4.0")}
    assert (u["repeat"], u["read_taxids"]) == ("yes", "1350:8,1578:1")
    assert {(m["contig"], m["repeat"]) for m in members} == {
        ("S1__contig_1", "no"), ("S2__contig_3", "yes"), ("S2__contig_5", "no")}


def test_cli_empty(tmp_path):
    cand, fasta = tmp_path / "c.tsv", tmp_path / "c.fa"
    cand.write_text(CAND_HEAD)
    fasta.write_text("")
    prefix = tmp_path / "u"
    assert pu.main(["--candidates", str(cand), "--fasta", str(fasta), "--isescan", "--recruit",
                    "--out-prefix", str(prefix)]) == 0
    assert open(str(prefix) + ".tsv").read() == "\t".join(pu.UNIT_COLUMNS) + "\n"
    assert open(str(prefix) + "_members.tsv").read() == "\t".join(pu.MEMBER_COLUMNS) + "\n"


def test_read_taxa_by_name_with_a_taxdump(tmp_path):
    a = seq(20000, 1)
    fasta = tmp_path / "c.fa"
    fasta.write_text(">S1__contig_1\n%s\n>S2__contig_3\n%s\n" % (a, a[3000:] + a[:3000]))
    cand = tmp_path / "c.tsv"
    cand.write_text(CAND_HEAD
                    + "S1__contig_1\tS1\t20000\t40\tno\t10\t0.9\t0.1\t1\t\trep27\t\tMOBV\t\tmobilizable\ttet(L)\t\t1350:6\n"
                    + "S2__contig_3\tS2\t20000\t20\tno\t5\t0.9\t0.1\t1\t\trep27\t\tMOBV\t\tmobilizable\ttet(L)\t\t1350:2,9999:1\n")
    tax = write_taxdump(tmp_path / "taxdump", TEST_TAXA)
    prefix = tmp_path / "u"
    assert pu.main(["--candidates", str(cand), "--fasta", str(fasta), "--taxdump", tax, "--out-prefix", str(prefix)]) == 0
    (u,) = list(csv.DictReader(open(str(prefix) + ".tsv"), delimiter="\t"))
    assert u["read_taxa"] == "Enterococcus:8,9999:1"          # a taxid missing from the taxdump stays a taxid
    assert u["read_taxids"] == "1350:8,9999:1"
    assert u["repeat"] == "no"
