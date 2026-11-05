"""bin_recovery.py: community bins with read-vote species, panel status and within-subject depth/breadth."""
import csv

from conftest import load_bin

br = load_bin("bin_recovery")
COV_HEAD = "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"


def w(path, text):
    path.write_text(text)
    return str(path)


def test_recover(tmp_path):
    tax = tmp_path / "taxdump"
    tax.mkdir()
    w(tax / "nodes.dmp", "".join("%s\t|\t%s\t|\t%s\t|\n" % r for r in [
        ("1", "1", "no rank"), ("724", "1", "genus"), ("727", "724", "species"), ("729", "724", "species"),
        ("99999", "727", "strain")]))
    w(tax / "names.dmp", "".join("%s\t|\t%s\t|\t\t|\tscientific name\t|\n" % r for r in [
        ("1", "root"), ("724", "Haemophilus"), ("727", "Haemophilus influenzae"),
        ("729", "Haemophilus parainfluenzae"), ("99999", "Haemophilus influenzae strain X")]))
    c2b = w(tmp_path / "contig2bin.tsv", "sample\tbin\tcontig\nS1\tS1.bin.1\tS1__c1\nS1\tS1.bin.1\tS1__c2\n"
            "S2\tS2.bin.1\tS2__c1\n")
    lengths = w(tmp_path / "len.tsv", "S1__c1\t1000\t\nS1__c2\t1000\t\nS2__c1\t2000\t\n")
    samples = w(tmp_path / "samples.tsv", "sample\tpatient\nS1\tP1\nS2\tP1\n")
    r2c = [w(tmp_path / "S1.read2contig.tsv", "r1\tS1__c1\nr2\tS1__c1\nr3\tS1__c2\nr4\tS2__c1\n"),
           w(tmp_path / "S2.read2contig.tsv", "q1\tS2__c1\n")]
    kr = [w(tmp_path / "S1.kraken2.classifiedreads.txt", "C\tr1\t727\t900\t-\nC\tr2\t99999\t900\t-\n"
            "C\tr3\t729\t900\t-\nU\tr5\t0\t900\t-\n"),
          w(tmp_path / "S2.kraken2.classifiedreads.txt", "C\tq1\tHaemophilus parainfluenzae (taxid 729)\t900\t-\n")]
    cov = [w(tmp_path / "S1.bins_coverage.tsv", COV_HEAD + "S1__c1\t1\t1000\t50\t1000\t100\t10\t20\t60\n"
             "S1__c2\t1\t1000\t5\t500\t50\t2\t20\t60\nS2__c1\t1\t2000\t2\t200\t10\t0.2\t20\t60\n"),
           w(tmp_path / "S2.bins_coverage.tsv", COV_HEAD + "S2__c1\t1\t2000\t200\t2000\t100\t30\t20\t60\n")]
    checkm2 = w(tmp_path / "quality_report.tsv", "Name\tCompleteness\tContamination\nS1.bin.1\t95.0\t2.0\n")
    cand = w(tmp_path / "cand.tsv", "sample\tspecies\tdecision\nS1\tHaemophilus influenzae\tbelow_est_depth\n")
    assert br.main(["--contig2bin", c2b, "--lengths", lengths, "--samples", samples, "--read2contig", *r2c,
                    "--kraken", *kr, "--coverage", *cov, "--checkm2", checkm2, "--candidates", cand,
                    "--taxdump", str(tax), "--outdir", str(tmp_path)]) == 0

    bins = {r["bin"]: r for r in csv.DictReader(open(tmp_path / "recovered_bins.tsv"), delimiter="\t")}
    b1, b2 = bins["S1.bin.1"], bins["S2.bin.1"]
    assert (b1["species"], b1["species_vote"], b1["genus"]) == ("Haemophilus influenzae", "0.667", "Haemophilus")
    assert (b1["size_bp"], b1["mimag_grade"], b1["panel_decision"]) == ("2000", "high", "below_est_depth")
    assert (b1["own_depth"], b1["own_breadth"], b1["n_reads"]) == ("6.00", "0.750", "3")
    assert (b2["species"], b2["panel_decision"], b2["mimag_grade"]) == ("Haemophilus parainfluenzae", "not_in_panel", "")
    assert list(b1) == br.BIN_COLUMNS

    per = {(r["bin"], r["read_sample"]): r for r in csv.DictReader(open(tmp_path / "bin_persistence.tsv"),
                                                                    delimiter="\t")}
    assert set(per) == {("S1.bin.1", "S1"), ("S2.bin.1", "S1"), ("S2.bin.1", "S2")}
    assert (per[("S2.bin.1", "S1")]["breadth"], per[("S2.bin.1", "S1")]["present"]) == ("0.100", "no")
    assert (per[("S2.bin.1", "S2")]["mean_depth"], per[("S2.bin.1", "S2")]["present"]) == ("30.00", "yes")


def test_subject_column_and_no_candidates(tmp_path):
    tax = tmp_path / "taxdump"
    tax.mkdir()
    w(tax / "nodes.dmp", "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tgenus\t|\n3\t|\t2\t|\tspecies\t|\n")
    w(tax / "names.dmp", "1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tLactobacillus\t|\t\t|\tscientific name\t|\n"
      "3\t|\tLactobacillus crispatus\t|\t\t|\tscientific name\t|\n")
    c2b = w(tmp_path / "contig2bin.tsv", "sample\tbin\tcontig\nA4\tA4.bin.1\tA4__c1\n")
    lengths = w(tmp_path / "len.tsv", "A4__c1\t1000\n")
    samples = w(tmp_path / "samples.tsv", "sample\tsubject\nA4\tPEN4\nA5\tPEN4\n")
    r2c = w(tmp_path / "A4.read2contig.tsv", "r1\tA4__c1\n")
    kr = w(tmp_path / "A4.kraken2.classifiedreads.txt", "C\tr1\t3\t900\t-\n")
    cov = [w(tmp_path / "A4.bins_coverage.tsv", COV_HEAD + "A4__c1\t1\t1000\t9\t1000\t100\t9\t20\t60\n"),
           w(tmp_path / "A5.bins_coverage.tsv", COV_HEAD + "A4__c1\t1\t1000\t3\t800\t80\t3\t20\t60\n")]
    assert br.main(["--contig2bin", c2b, "--lengths", lengths, "--samples", samples, "--read2contig", r2c,
                    "--kraken", kr, "--coverage", *cov, "--taxdump", str(tax), "--outdir", str(tmp_path)]) == 0
    b = list(csv.DictReader(open(tmp_path / "recovered_bins.tsv"), delimiter="\t"))[0]
    assert (b["subject"], b["species"], b["panel_decision"]) == ("PEN4", "Lactobacillus crispatus", "")
    per = list(csv.DictReader(open(tmp_path / "bin_persistence.tsv"), delimiter="\t"))
    assert [(r["subject"], r["read_sample"], r["present"]) for r in per] == [("PEN4", "A4", "yes"), ("PEN4", "A5", "yes")]
