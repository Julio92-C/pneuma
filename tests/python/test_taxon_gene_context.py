"""taxon_gene_context.py: every read hit of every context database with taxon fields and a contig host vote."""
import csv

import pytest

from conftest import TEST_TAXA, load_bin, write_taxdump

ec = load_bin("extend_context")
tg = load_bin("taxon_gene_context")
br = load_bin("bin_recovery")


def hc_row(**kw):
    row = dict.fromkeys(ec.COLUMNS, "")
    row.update(sample="S1", database="card", gene="tet(M)", contig_context="no_contig_match", final_context="unresolved",
               context_source="read", context_conflict="no", intrinsic_variant="no")
    row.update(kw)
    return row


def write_rows(path, columns, rows):
    with open(str(path), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return str(path)


def read(path):
    return list(csv.DictReader(open(str(path)), delimiter="\t"))


def test_lineage_levels(tmp_path):
    tax = br.read_taxonomy(write_taxdump(tmp_path / "t", TEST_TAXA))
    ec_ = tg.lineage("562", tax)
    assert (ec_["kraken_rank"], ec_["kraken_name"], ec_["species"], ec_["genus"], ec_["family"], ec_["phylum"]) == (
        "species", "Escherichia coli", "Escherichia coli", "Escherichia", "Enterobacteriaceae", "Pseudomonadota")
    assert ec_["taxon_level"] == "species"
    assert tg.lineage("1578", tax)["taxon_level"] == "genus"
    assert tg.lineage("1", tax)["taxon_level"] == "above_family"
    assert tg.lineage("0", tax)["taxon_level"] == "unclassified"
    assert tg.lineage("", tax)["taxon_level"] == "unclassified"
    unknown = tg.lineage("999999", tax)
    assert unknown["taxon_level"] == "unclassified"
    assert (unknown["kraken_name"], unknown["species"], unknown["genus"], unknown["family"], unknown["phylum"]) == (
        "", "", "", "", "")


def test_host_vote_and_weak_flag(tmp_path):
    hc = write_rows(tmp_path / "hit_context_card.tsv", ec.COLUMNS, [
        hc_row(seq="S1__r1", contig="S1__contig_1", contig_context="plasmid", final_context="plasmid"),
        hc_row(seq="S1__r1", gene="sul1", contig="S1__contig_1", contig_context="plasmid", final_context="plasmid"),
        hc_row(seq="S1__r2", contig="S1__contig_1", contig_context="plasmid", final_context="plasmid"),
        hc_row(seq="S1__r3", contig="S1__contig_1", contig_context="plasmid", final_context="plasmid"),
        hc_row(seq="S1__r4", contig="S1__contig_2", contig_context="chromosome", final_context="chromosome"),
        hc_row(seq="S1__r5")])
    taxa = tmp_path / "taxa.tsv"
    taxa.write_text("S1__r1\t562\nS1__r2\t562\nS1__r3\t1351\nS1__r4\t562\nS1__r5\t562\n")
    out, res = tmp_path / "hit_context.tsv", tmp_path / "res.tsv"
    assert tg.main(["--hit-context", hc, "--taxa", str(taxa), "--taxdump", write_taxdump(tmp_path / "t", TEST_TAXA),
                    "--out", str(out), "--resolution", str(res)]) == 0
    assert open(str(out)).readline().rstrip("\n").split("\t") == tg.COLUMNS
    rows = read(out)
    assert len(rows) == 6
    r1 = rows[0]
    assert (r1["element"], r1["species"], r1["contig_taxon"]) == ("AMR", "Escherichia coli", "Escherichia coli")
    # contig_1: three distinct hit reads (r1 counted once), two E. coli -> 0.667, not weak
    assert (r1["contig_vote"], r1["contig_vote_reads"], r1["weak_vote"]) == ("0.667", "2", "no")
    r4 = next(r for r in rows if r["seq"] == "S1__r4")
    assert (r4["contig_vote"], r4["contig_vote_reads"], r4["weak_vote"]) == ("1.000", "1", "yes")   # one read only
    r5 = next(r for r in rows if r["seq"] == "S1__r5")
    assert (r5["contig_taxon"], r5["contig_vote"], r5["weak_vote"]) == ("", "", "NA")
    all_amr = next(r for r in read(res) if (r["element"], r["context"]) == ("AMR", "all"))
    assert (all_amr["hits"], all_amr["species_pct"]) == ("6", "100.0")


def test_a_tied_host_vote_is_weak():
    rows = [dict(seq="S1__r%d" % i, contig="S1__contig_1") for i in range(1, 5)]
    fields = {"S1__r1": {"species": "Escherichia coli", "genus": "Escherichia"},
              "S1__r2": {"species": "Escherichia coli", "genus": "Escherichia"},
              "S1__r3": {"species": "Klebsiella pneumoniae", "genus": "Klebsiella"},
              "S1__r4": {"species": "Klebsiella pneumoniae", "genus": "Klebsiella"}}
    # 2 vs 2: share 0.5 and 2 reads pass the share and read floors, but the tie makes the vote weak
    assert tg.host_votes(rows, fields)["S1__contig_1"] == ("Escherichia coli", 0.5, 2, True)
    fields["S1__r4"] = {"species": "Escherichia coli", "genus": "Escherichia"}
    assert tg.host_votes(rows, fields)["S1__contig_1"] == ("Escherichia coli", 0.75, 3, False)


def test_read_loci_fill_in_a_database_without_a_hit_context_table(tmp_path):
    head = "sample\troute\tseq\tseq_len\tgene\tdrug_class\tstart\tend\tcategory\treason\n"
    hc = write_rows(tmp_path / "hit_context_card.tsv", ec.COLUMNS, [hc_row(seq="S1__r1")])
    card = tmp_path / "loci_reads_card.tsv"
    card.write_text(head + "S1\treads\tS1__r1\t9000\ttet(M)\tx\t1\t900\tunresolved\tshort_flank\n")
    vfdb = tmp_path / "loci_reads_vfdb.tsv"
    vfdb.write_text(head + "S1\treads\tS1__r9\t9000\tfimH\t\t1\t900\tchromosome\tgenomad\n")
    out, res = tmp_path / "o.tsv", tmp_path / "r.tsv"
    assert tg.main(["--hit-context", hc, "--read-loci", str(card), str(vfdb),
                    "--taxdump", write_taxdump(tmp_path / "t", TEST_TAXA), "--out", str(out), "--resolution", str(res)]) == 0
    rows = read(out)
    assert [(r["database"], r["seq"]) for r in rows] == [("card", "S1__r1"), ("vfdb", "S1__r9")]   # card not doubled
    vf = rows[1]
    assert (vf["element"], vf["contig_context"], vf["final_context"], vf["taxon_level"], vf["weak_vote"]) == (
        "VFs", "no_contig_match", "chromosome", "unclassified", "NA")


def test_score_only_tables_are_refused(tmp_path):
    f = tmp_path / "loci_reads_card_score_only.tsv"
    f.write_text("sample\troute\tseq\tseq_len\tgene\tdrug_class\tstart\tend\tcategory\treason\n")
    with pytest.raises(ValueError):
        tg.load([], [str(f)])
