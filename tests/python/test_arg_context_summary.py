"""arg_context_summary.py: per-sample / per-class / per-group shares of ARG contexts, and route concordance."""
import csv

from conftest import load_bin

sm = load_bin("arg_context_summary")
HEAD = "\t".join(["sample", "route", "seq", "seq_len", "gene", "drug_class", "start", "end", "category"]) + "\n"


def loci_file(tmp_path, name, rows):
    p = tmp_path / name
    p.write_text(HEAD + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(p)


def read(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def test_summary_tables(tmp_path):
    reads = loci_file(tmp_path, "reads.tsv", [
        ("S1", "reads", "S1__r1", 9000, "tet(M)", "tetracycline", 3000, 4000, "chromosomal_mge"),
        ("S1", "reads", "S1__r2", 9000, "tet(M)", "tetracycline", 3000, 4000, "chromosomal_mge"),
        ("S1", "reads", "S1__r3", 9000, "blaTEM-1", "beta-lactam;penam", 3000, 4000, "plasmid"),
        ("S2", "reads", "S2__r1", 9000, "mef(A)", "macrolide", 3000, 4000, "unresolved"),
    ])
    contigs = loci_file(tmp_path, "contigs.tsv", [
        ("S1", "contigs", "S1__contig_1", 50000, "tet(M)", "tetracycline", 101, 200, "chromosome"),
        ("S1", "contigs", "S1__contig_2", 6000, "blaTEM-1", "beta-lactam;penam", 1, 100, "plasmid"),
    ])
    bedcov = tmp_path / "S1.bedcov"
    bedcov.write_text("S1__contig_1\t100\t200\ttet(M)\t300\nS1__contig_2\t0\t100\tblaTEM-1\t900\n")
    meta = tmp_path / "meta.csv"
    meta.write_text("sample,sample_type\nS1,Sputum\nS2,Exhale\n")
    assert sm.main(["--loci", reads, contigs, "--bedcov", str(bedcov), "--metadata", str(meta),
                    "--group-cols", "sample_type", "--outdir", str(tmp_path)]) == 0

    s = {(r["sample"], r["route"], r["category"]): r for r in read(tmp_path / "arg_context_by_sample.tsv")}
    assert s[("S1", "reads", "chromosomal_mge")]["fraction"] == "0.6667"
    assert s[("S1", "contigs", "plasmid")]["weight"] == "9.00"       # 900 depth-bp over 100 bp
    assert s[("S1", "contigs", "plasmid")]["fraction"] == "0.7500"   # 9 / (9 + 3)
    assert s[("S2", "reads", "unresolved")]["fraction"] == "1.0000"

    c = {(r["route"], r["drug_class"], r["category"]): r for r in read(tmp_path / "arg_context_by_class.tsv")}
    assert c[("reads", "penam", "plasmid")]["n_loci"] == "1"

    k = {(r["sample"], r["gene"]): r for r in read(tmp_path / "arg_context_concordance.tsv")}
    assert (k[("S1", "blaTEM-1")]["agree"], k[("S1", "tet(M)")]["agree"]) == ("yes", "no")
    assert k[("S2", "mef(A)")]["agree"] == "NA"

    g = {(r["group"], r["route"], r["category"]): r for r in read(tmp_path / "arg_context_by_group.tsv")}
    assert g[("Sputum", "reads", "plasmid")]["pooled_fraction"] == "0.3333"
    assert g[("Exhale", "reads", "unresolved")]["n_samples"] == "1"
