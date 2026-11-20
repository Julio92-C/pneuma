"""arg_context_summary.py: shares of the final genomic context per sample, drug class, metadata group and gene."""
import csv

from conftest import load_bin

sm = load_bin("arg_context_summary")
VARIANT = "Bifidobacterium_adolescentis_rpoB_mutants_conferring_resistance_to_rifampicin"
COLS = ["sample", "seq", "database", "gene", "drug_class", "final_context", "intrinsic_variant"]
HITS = [("S1", "S1__r1", "card", "tet(M)", "tetracycline", "plasmid", "no"),
        ("S1", "S1__r2", "card", "tet(M)", "tetracycline", "chromosome", "no"),
        ("S1", "S1__r3", "card", VARIANT, "rifamycin", "chromosome", "yes"),
        ("S1", "S1__r4", "card", "blaTEM-1", "beta-lactam;penam", "plasmid", "no"),
        ("S2", "S2__r1", "card", "mef(A)", "macrolide", "unresolved", "no"),
        ("S2", "S2__r2", "vfdb", "fimH", "", "chromosome", "no"),
        ("S3", "S3__r1", "card", "tet(M)", "tetracycline", "plasmid", "no")]


def hit_context(tmp_path):
    p = tmp_path / "hit_context.tsv"
    p.write_text("\t".join(COLS) + "\n" + "".join("\t".join(h) + "\n" for h in HITS))
    return str(p)


def read(path):
    return list(csv.DictReader(open(str(path)), delimiter="\t"))


def test_by_sample_class_and_gene(tmp_path):
    assert sm.main(["--hit-context", hit_context(tmp_path), "--outdir", str(tmp_path / "s")]) == 0
    s = {(r["sample"], r["database"], r["category"]): r for r in read(tmp_path / "s" / "arg_context_by_sample.tsv")}
    assert (s[("S1", "card", "plasmid")]["n_hits"], s[("S1", "card", "plasmid")]["fraction"]) == ("2", "0.5000")
    chrom = s[("S1", "card", "chromosome")]
    assert (chrom["n_hits"], chrom["n_intrinsic"], chrom["intrinsic_share"]) == ("2", "1", "0.5000")
    assert (s[("S1", "card", "chromosomal_mge")]["fraction"], s[("S1", "card", "chromosomal_mge")]["intrinsic_share"]) == (
        "0.0000", "NA")
    assert s[("S2", "vfdb", "chromosome")]["fraction"] == "1.0000"

    c = {(r["drug_class"], r["category"]): r for r in read(tmp_path / "s" / "arg_context_by_class.tsv")}
    assert c[("tetracycline", "plasmid")]["fraction"] == "0.6667"
    assert c[("penam", "plasmid")]["n_hits"] == "1"
    assert c[("rifamycin", "chromosome")]["n_intrinsic"] == "1"
    assert not any(k[0] in ("", "unknown") for k in c)                      # VFDB hits have no drug class row

    g = {(r["database"], r["gene"]): r for r in read(tmp_path / "s" / "gene_context.tsv")}
    tet = g[("card", "tet(M)")]
    assert (tet["n_hits"], tet["n_samples"], tet["plasmid"], tet["chromosome"], tet["n_intrinsic"]) == ("3", "2", "2", "1", "0")
    assert not (tmp_path / "s" / "arg_context_by_group.tsv").exists()


def test_by_group(tmp_path):
    meta = tmp_path / "meta.csv"
    meta.write_text("sample,sample_type\nS1,Sputum\nS2,Exhale\nS3,Sputum\n")
    assert sm.main(["--hit-context", hit_context(tmp_path), "--metadata", str(meta), "--group-cols", "sample_type",
                    "--outdir", str(tmp_path)]) == 0
    g = {(r["group"], r["database"], r["category"]): r for r in read(tmp_path / "arg_context_by_group.tsv")}
    sp = g[("Sputum", "card", "plasmid")]
    assert (sp["n_samples"], sp["n_hits"], sp["pooled_fraction"], sp["mean_sample_fraction"]) == ("2", "3", "0.6000", "0.7500")
    assert (g[("Exhale", "card", "unresolved")]["n_samples"], g[("Exhale", "card", "unresolved")]["pooled_fraction"]) == (
        "1", "1.0000")
    assert g[("Exhale", "vfdb", "chromosome")]["n_hits"] == "1"


def test_missing_group_column_fails(tmp_path, capsys):
    meta = tmp_path / "meta.csv"
    meta.write_text("sample,sample_type\nS1,Sputum\n")
    assert sm.main(["--hit-context", hit_context(tmp_path), "--metadata", str(meta), "--group-cols", "ward",
                    "--outdir", str(tmp_path)]) == 1
    assert "ward" in capsys.readouterr().err
