"""is_flanks.py: IS elements next to ARGs on contigs; composite transposons."""
import csv

from conftest import load_bin

isf = load_bin("is_flanks")

ABR = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY"
       "\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n")
ISE = "seqID\tfamily\tcluster\tisBegin\tisEnd\tisLen\tstrand\ttype\n"


def write(path, text):
    path.write_text(text)
    return str(path)


def hit(seq, start, end, gene, db="card"):
    return "c.fa\t%s\t%d\t%d\t+\t%s\t1-1/1\t===\t0/0\t100\t99\t%s\tX\tp\tr\n" % (seq, start, end, gene, db)


def test_flanks_and_composite(tmp_path):
    hits = write(tmp_path / "h.tsv", ABR + hit("S1__contig_1", 5000, 6000, "aph(3')-Ia") + hit("S1__contig_2", 5000, 6000, "tet(M)")
                 + hit("S1__contig_3", 5000, 6000, "ErmB") + hit("S1__contig_1", 100, 400, "repA", "plasmidfinder"))
    ise = write(tmp_path / "i.tsv", ISE
                + "S1__contig_1\tIS6\tIS26\t3800\t4800\t1001\t+\tc\nS1__contig_1\tIS6\tIS26\t6300\t7100\t801\t-\tc\n"
                + "S1__contig_2\tIS6\tIS26\t3800\t4800\t1001\t+\tc\nS1__contig_2\tIS6\tIS26\t6300\t7100\t801\t+\tc\n"
                + "S1__contig_3\tIS3\tIS3\t1000\t2000\t1001\t+\tc\n")
    out = tmp_path / "o.tsv"
    assert isf.main(["--isescan", ise, "--hits", hits, "--window", "500", "--out", str(out)]) == 0
    rows = {r["contig"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert len(rows) == 3                                         # card hits only, one row each
    c1 = rows["S1__contig_1"]
    assert (c1["sample"], c1["left_family"], c1["left_distance"], c1["right_family"], c1["right_distance"]) == (
        "S1", "IS6", "200", "IS6", "300")
    assert c1["composite_transposon"] == "yes"
    assert c1["orientation"] == "inverted"                        # inverted: different strands (+ vs -)
    assert rows["S1__contig_2"]["composite_transposon"] == "yes"  # same family, now allowed in any orientation
    assert rows["S1__contig_2"]["orientation"] == "direct"        # direct: same strands (+ vs +)
    c3 = rows["S1__contig_3"]
    assert (c3["left_is"], c3["right_is"], c3["composite_transposon"]) == ("", "", "no")   # 3 kb away
    assert c3["orientation"] == ""                                # empty when not composite


def test_composite_with_unknown_strand(tmp_path):
    """Composite with same family but one strand missing -> orientation='unknown'."""
    hits = write(tmp_path / "h.tsv", ABR + hit("S1__contig_1", 5000, 6000, "aph(3')-Ia"))
    ise = write(tmp_path / "i.tsv", ISE
                + "S1__contig_1\tIS6\tIS26\t3800\t4800\t1001\t+\tc\nS1__contig_1\tIS6\tIS26\t6300\t7100\t801\t\tc\n")
    out = tmp_path / "o.tsv"
    assert isf.main(["--isescan", ise, "--hits", hits, "--window", "500", "--out", str(out)]) == 0
    rows = {r["contig"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    c1 = rows["S1__contig_1"]
    assert c1["composite_transposon"] == "yes"
    assert c1["orientation"] == "unknown"


def test_empty_inputs(tmp_path):
    out = tmp_path / "o.tsv"
    assert isf.main(["--isescan", write(tmp_path / "i.tsv", ""), "--hits", write(tmp_path / "h.tsv", ABR),
                     "--out", str(out)]) == 0
    assert out.read_text() == "\t".join(isf.COLUMNS) + "\n"
