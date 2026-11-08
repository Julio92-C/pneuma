"""pathogen_candidates.py: one decision per sample x panel species."""
import csv

import pytest

from conftest import DATA, load_bin

pc = load_bin("pathogen_candidates")
MINI = DATA / "pathogen_mini"
PANEL_HEADER = "\t".join(load_bin("resolve_panel").OUT_COLUMNS)


def write(path, text):
    path.write_text(text)
    return path


def panel(tmp_path):
    rows = [
        ["Streptococcus pneumoniae", "Streptococcus_pneumoniae", "1", "", "reconstruct", "1313", "Streptococcus pneumoniae",
         "scientific name", "1313;170187", "GCF_1", "", "", "95.0", "Streptococcus_pneumoniae"],
        ["Mycobacterium tuberculosis", "Mycobacterium_tuberculosis", "D", "", "detect", "1773", "Mycobacterium tuberculosis",
         "scientific name", "1773", "", "", "", "", ""],
        ["Haemophilus influenzae", "Haemophilus_influenzae", "1", "", "reconstruct", "727", "Haemophilus influenzae",
         "scientific name", "727", "GCF_2", "", "", "95.0", ""],
    ]
    return write(tmp_path / "panel_resolved.tsv", PANEL_HEADER + "\n" + "".join("\t".join(r) + "\n" for r in rows))


def kreport(tmp_path, sample, spn, mtb, unclassified=10):
    lines = ["%6.2f\t%d\t%d\tU\t0\tunclassified" % (0, unclassified, unclassified),
             "%6.2f\t%d\t0\tR\t1\troot" % (0, spn + mtb),
             "%6.2f\t%d\t%d\tS\t1313\t    Streptococcus pneumoniae" % (0, spn, spn),
             "%6.2f\t%d\t%d\tS\t1773\t    Mycobacterium tuberculosis" % (0, mtb, mtb)]
    return write(tmp_path / ("%s.kraken2.report.txt" % sample), "\n".join(lines) + "\n")


def test_decisions(tmp_path):
    samples = write(tmp_path / "samples.tsv", "sample\tpatient\tis_control\tlow_depth\tn_reads\n"
                    "S1\tP1\tfalse\tfalse\t1000\nS2\t\tfalse\ttrue\t40\nS3\tP1\tfalse\tfalse\t500\n")
    reports = [kreport(tmp_path, "S1", 120, 3), kreport(tmp_path, "S2", 30, 0), kreport(tmp_path, "S3", 4, 0)]
    bracken = [write(tmp_path / "S1.tsv", "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
                     "new_est_reads\tfraction_total_reads\nStreptococcus pneumoniae\t1313\tS\t120\t5\t125\t0.97\n")]
    out, sout = tmp_path / "c.tsv", tmp_path / "s.tsv"
    rc = pc.main(["--panel", str(panel(tmp_path)), "--samples", str(samples), "--kreports", *map(str, reports),
                  "--bracken", *map(str, bracken), "--min-reads", "10", "--out", str(out), "--samples-out", str(sout)])
    assert rc == 0
    rows = {(r["sample"], r["species"]): r for r in csv.DictReader(open(out), delimiter="\t")}
    assert len(rows) == 9
    s1 = rows[("S1", "Streptococcus pneumoniae")]
    assert (s1["decision"], s1["clade_reads"], s1["bracken_reads"]) == ("extract", "120", "125")
    assert rows[("S1", "Mycobacterium tuberculosis")]["decision"] == "detect"
    assert rows[("S1", "Haemophilus influenzae")]["decision"] == "not_detected"
    assert rows[("S1", "Haemophilus influenzae")]["bracken_reads"] == "0"
    assert all(rows[("S2", sp)]["decision"] == "low_depth" for sp in
               ("Streptococcus pneumoniae", "Mycobacterium tuberculosis", "Haemophilus influenzae"))
    assert rows[("S2", "Streptococcus pneumoniae")]["bracken_reads"] == "NA"
    assert rows[("S3", "Streptococcus pneumoniae")]["decision"] == "below_min_reads"
    assert open(out).readline().rstrip("\n").split("\t") == pc.CAND_COLUMNS
    srows = {r["sample"]: r for r in csv.DictReader(open(sout), delimiter="\t")}
    assert srows["S1"]["unclassified_pct"] == "7.52" and srows["S1"]["classified_reads"] == "123"


def test_missing_kreport_fails(tmp_path):
    samples = write(tmp_path / "samples.tsv", "sample\tpatient\tis_control\tlow_depth\tn_reads\nS9\t\tfalse\tfalse\t5\n")
    with pytest.raises(ValueError, match="no Kraken2 report for sample S9"):
        pc.run(panel(tmp_path), samples, [], [], 10, tmp_path / "c.tsv", tmp_path / "s.tsv")


def test_real_mini_fixture(tmp_path):
    """Contract: the Kraken2/Bracken files written like PNEUMA's TAXONOMY outputs."""
    rp = load_bin("resolve_panel")
    resolved = tmp_path / "panel_resolved.tsv"
    rp.run(MINI / "panel.csv", MINI / "panel_references.tsv", DATA / "pathogen_fixture" / "taxdump", resolved)
    out = tmp_path / "c.tsv"
    pc.run(resolved, MINI / "samples.tsv", sorted(MINI.glob("*.kraken2.report.txt")),
           sorted((MINI / "bracken").glob("*.tsv")), 10, out, tmp_path / "s.tsv")
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [r["decision"] for r in rows] == ["extract", "extract", "extract"]
    assert int(rows[0]["clade_reads"]) > int(rows[1]["clade_reads"]) > 0


# --- estimated-depth filter (--min-est-depth) -------------------------------------------------
import gzip


def refdir(tmp_path, length=10_000):
    d = tmp_path / "refs"
    d.mkdir()
    with gzip.open(d / "GCF_1.fna.gz", "wt") as fh:      # S. pneumoniae reference
        fh.write(">c1\n" + "A" * (length // 2) + "\n>c2\n" + "C" * (length // 2) + "\n")
    return d


def stats(tmp_path, sample, n, bases):
    d = tmp_path / "seqkit"
    d.mkdir(exist_ok=True)
    return write(d / ("%s.tsv" % sample), "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\n"
                 "%s.fq.gz\tFASTQ\tDNA\t%d\t%d\t1\t1\t1\n" % (sample, n, bases))


def run_depth(tmp_path, rows, stat_rows, min_est_depth=2.5, with_refs=True, gate_controls=False):
    samples = write(tmp_path / "samples.tsv", "sample\tpatient\tis_control\tlow_depth\tn_reads\n" +
                    "".join("%s\t%s\t%s\tfalse\t1000\n" % (s, p, c) for s, p, c, _ in rows))
    reports = [kreport(tmp_path, s, spn, 0) for s, _, _, spn in rows]
    st = [stats(tmp_path, s, n, b) for s, n, b in stat_rows]
    args = ["--panel", str(panel(tmp_path)), "--samples", str(samples), "--kreports", *map(str, reports),
            "--min-reads", "10", "--min-est-depth", str(min_est_depth),
            "--out", str(tmp_path / "c.tsv"), "--samples-out", str(tmp_path / "s.tsv")]
    if st:
        args += ["--seqkit-stats", *map(str, st)]
    if with_refs:
        args += ["--refdir", str(refdir(tmp_path))]
    if gate_controls:
        args += ["--gate-controls"]
    assert pc.main(args) == 0
    out = csv.DictReader(open(tmp_path / "c.tsv"), delimiter="\t")
    return {r["sample"]: r for r in out if r["species"] == "Streptococcus pneumoniae"}


def test_est_depth_filters_hopeless_candidates(tmp_path):
    # mean read length 1000 bp, genome 10 kb: depth = clade reads / 10
    rows = run_depth(tmp_path,
                     [("A", "P1", "false", 120),    # 12x  -> extract
                      ("B", "P2", "false", 20),     #  2x, patient P2 total 4x -> extract (pooling)
                      ("C", "P2", "false", 20),
                      ("D", "", "false", 20),       #  2x, no patient -> below_est_depth
                      ("E", "P3", "true", 20)],     #  control: not extracted (no --gate-controls)
                     [(s, 1000, 1_000_000) for s in "ABCDE"])
    assert rows["A"]["decision"] == "extract"
    assert rows["B"]["decision"] == rows["C"]["decision"] == "extract"
    assert rows["D"]["decision"] == "below_est_depth"
    assert rows["E"]["decision"] == "control"
    assert float(rows["D"]["est_depth"]) == pytest.approx(2.0)


def test_est_depth_unknown_keeps_old_behaviour(tmp_path):
    rows = run_depth(tmp_path, [("D", "", "false", 20)], [])                  # no seqkit stats
    assert rows["D"]["decision"] == "extract" and rows["D"]["est_depth"] == "NA"


def test_est_depth_zero_disables_filter(tmp_path):
    rows = run_depth(tmp_path, [("D", "", "false", 20)], [("D", 1000, 1_000_000)], min_est_depth=0)
    assert rows["D"]["decision"] == "extract"


def test_control_is_not_extracted_by_default(tmp_path):
    rows = run_depth(tmp_path, [("K", "", "true", 500)], [("K", 1000, 1_000_000)])
    assert rows["K"]["decision"] == "control"


def test_control_extracted_with_gate_controls(tmp_path):
    rows = run_depth(tmp_path, [("K", "", "true", 500)], [("K", 1000, 1_000_000)], gate_controls=True)
    assert rows["K"]["decision"] == "extract"


def test_genome_size_column(tmp_path):
    rows = run_depth(tmp_path, [("A", "P1", "false", 120)], [("A", 1000, 1_000_000)])
    assert rows["A"]["genome_size"] == "10000"
    rows = run_depth(tmp_path / "noref", [("A", "P1", "false", 120)], [], with_refs=False) if (tmp_path / "noref").mkdir() is None else None
    assert rows["A"]["genome_size"] == ""
