"""arg_context.py: genomic context (plasmid / chromosomal_mge / chromosome / unresolved) of ARG hits."""
import csv

import pytest
from conftest import load_bin

ac = load_bin("arg_context")

ABR_HEAD = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY"
            "\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n")


def abr_row(seq, start, end, gene, db="card", res="macrolide"):
    return "f.fa\t%s\t%d\t%d\t+\t%s\t1-1/1\t===\t0/0\t100.00\t99.00\t%s\tX1\tp\t%s\n" % (seq, start, end, gene, db, res)


def write(path, text):
    path.write_text(text)
    return str(path)


@pytest.mark.parametrize("kw,expected", [
    (dict(scores=(0.9, 0.05, 0.0)), ("plasmid", "genomad")),
    (dict(scores=(0.1, 0.9, 0.0), mob_molecule="plasmid"), ("unresolved", "conflict_mob_plasmid_genomad_chromosome")),
    (dict(scores=(0.4, 0.5, 0.0), mob_molecule="plasmid"), ("plasmid", "mobsuite")),
    (dict(scores=(0.9, 0.0, 0.0), mob_molecule="plasmid"), ("plasmid", "genomad+mobsuite")),
    (dict(scores=(0.0, 0.95, 0.0), mob_molecule="chromosome"), ("chromosome", "genomad")),
    (dict(scores=(0.0, 0.95, 0.0), mge_near=True), ("chromosomal_mge", "mge_genes_nearby")),
    (dict(scores=(0.1, 0.2, 0.7), in_provirus=True), ("chromosomal_mge", "provirus")),
    (dict(scores=(0.3, 0.4, 0.3)), ("unresolved", "low_scores")),
    (dict(scores=None, mob_molecule="chromosome"), ("unresolved", "no_scores")),
    (dict(scores=(0.9, 0.0, 0.0), flank_ok=False), ("unresolved", "short_flank")),
    (dict(scores=(0.2, 0.3, 0.0), mge_near=True), ("unresolved", "low_scores")),
])
def test_classify(kw, expected):
    assert ac.classify(**kw) == expected


def test_read_abricate_keeps_card_and_handles_repeated_headers_and_csv(tmp_path):
    tab = write(tmp_path / "a.tsv", ABR_HEAD + abr_row("S1__r1", 900, 100, "mef(A)") + ABR_HEAD +
                abr_row("S1__r2", 5, 50, "vfX", db="vfdb"))
    assert ac.read_abricate(tab) == [{"seq": "S1__r1", "start": 100, "end": 900, "gene": "mef(A)",
                                      "drug_class": "macrolide"}]
    csvf = write(tmp_path / "a.csv", ABR_HEAD.replace("\t", ",") + abr_row("S2__c1", 10, 20, "tet(M)").replace("\t", ","))
    assert ac.read_abricate(csvf)[0]["seq"] == "S2__c1"


def test_unprefixed_ids_fail():
    with pytest.raises(ValueError, match="sample"):
        ac.sample_of("read_without_prefix")


def test_main_reads_route(tmp_path):
    abr = write(tmp_path / "abr.tsv", ABR_HEAD +
                abr_row("S1__plas", 3000, 4000, "blaTEM-1", res="beta-lactam") +
                abr_row("S1__ice", 3000, 4000, "tet(M)", res="tetracycline") +
                abr_row("S2__chr", 3000, 4000, "mef(A)") +
                abr_row("S2__short", 200, 1200, "msr(D)"))
    lengths = write(tmp_path / "len.tsv", "S1__plas\t8000\t\nS1__ice\t9000\t\nS2__chr\t9000\t\nS2__short\t5000\t\n")
    scores = write(tmp_path / "agg.tsv", "seq_name\tchromosome_score\tplasmid_score\tvirus_score\n"
                   "S1__plas\t0.05\t0.93\t0.02\nS1__ice\t0.88\t0.10\t0.02\nS2__chr\t0.91\t0.05\t0.04\n"
                   "S2__short\t0.10\t0.85\t0.05\n")
    virus = write(tmp_path / "vir.tsv", "seq_name\tlength\ttopology\nS2__chr|provirus_100_900\t800\tProvirus\n")
    genes = write(tmp_path / "genes.tsv", "gene\tstart\tend\tannotation_conjscan\tannotation_description\n"
                  "S1__ice_1\t6000\t7000\tNA\tsite-specific integrase\n"
                  "S2__chr_1\t3000\t4000\tNA\tMFS efflux\n")
    mobt = write(tmp_path / "mobt.tsv", "sample_id\trep_type(s)\trelaxase_type(s)\nS1__plas\trep_cluster_1\t-\n"
                 "S1__ice\t-\tMOBV\n")
    out = tmp_path / "loci.tsv"
    assert ac.main(["--route", "reads", "--abricate", abr, "--lengths", lengths, "--genomad-scores", scores,
                    "--genomad-virus", virus, "--genomad-genes", genes, "--mob-typer", mobt,
                    "--out", str(out)]) == 0
    rows = {r["seq"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert (rows["S1__plas"]["category"], rows["S1__plas"]["reason"]) == ("plasmid", "genomad+mobsuite")
    assert (rows["S1__ice"]["category"], rows["S1__ice"]["mge_genes"]) == ("chromosomal_mge", "site-specific integrase")
    assert (rows["S2__chr"]["category"], rows["S2__chr"]["in_provirus"]) == ("chromosome", "no")
    assert (rows["S2__short"]["category"], rows["S2__short"]["reason"]) == ("unresolved", "short_flank")
    assert rows["S2__chr"]["sample"] == "S2" and list(rows["S2__chr"]) == ac.LOCUS_COLUMNS


def test_main_contigs_route_uses_mob_recon_and_proviruses(tmp_path):
    abr = write(tmp_path / "abr.tsv", ABR_HEAD + abr_row("S1__contig_1", 50, 900, "aph(6)-Id", res="aminoglycoside") +
                abr_row("S1__contig_2", 20000, 21000, "erm(B)"))
    lengths = write(tmp_path / "len.tsv", "S1__contig_1\t4000\nS1__contig_2\t300000\n")
    scores = write(tmp_path / "agg.tsv", "seq_name\tchromosome_score\tplasmid_score\tvirus_score\n"
                   "S1__contig_1\t0.30\t0.55\t0.15\nS1__contig_2\t0.60\t0.05\t0.35\n")
    virus = write(tmp_path / "vir.tsv", "seq_name\tlength\nS1__contig_2|provirus_15000_40000\t25001\n")
    mob = write(tmp_path / "mob.tsv", "sample_id\tmolecule_type\tprimary_cluster_id\tcontig_id\n"
                "S1\tplasmid\tAA123\tS1__contig_1 len=4000\nS1\tchromosome\t-\tS1__contig_2\n"
                "sample_id\tmolecule_type\tprimary_cluster_id\tcontig_id\n")
    out = tmp_path / "loci.tsv"
    assert ac.main(["--route", "contigs", "--abricate", abr, "--lengths", lengths, "--genomad-scores", scores,
                    "--genomad-virus", virus, "--mob-contigs", mob, "--min-flank", "0", "--out", str(out)]) == 0
    rows = {r["seq"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert (rows["S1__contig_1"]["category"], rows["S1__contig_1"]["mob_detail"]) == ("plasmid", "AA123")
    assert (rows["S1__contig_2"]["category"], rows["S1__contig_2"]["reason"]) == ("chromosomal_mge", "provirus")


def test_abricate_db_option_selects_vfdb_hits(tmp_path):
    abr = write(tmp_path / "abr.tsv", ABR_HEAD + abr_row("S1__r1", 3000, 4000, "tet(M)") +
                abr_row("S1__r1", 5000, 6000, "fimH", db="vfdb", res=""))
    lengths = write(tmp_path / "len.tsv", "S1__r1\t9000\n")
    scores = write(tmp_path / "agg.tsv", "seq_name\tchromosome_score\tplasmid_score\tvirus_score\n"
                   "S1__r1\t0.9\t0.05\t0.05\n")
    out = tmp_path / "loci.tsv"
    assert ac.main(["--route", "reads", "--abricate", abr, "--abricate-db", "VFDB", "--lengths", lengths,
                    "--genomad-scores", scores, "--out", str(out)]) == 0
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [(r["gene"], r["category"]) for r in rows] == [("fimH", "chromosome")]


@pytest.mark.parametrize("kw,expected", [
    (dict(scores=(0.9, 0.0, 0.0), plasmid_supported=False), ("unresolved", "plasmid_score_only")),
    (dict(scores=(0.9, 0.0, 0.0), plasmid_supported=False, mob_molecule="plasmid"), ("plasmid", "genomad+mobsuite")),
    (dict(scores=(0.9, 0.0, 0.0), plasmid_supported=True), ("plasmid", "genomad")),
    (dict(scores=(0.9, 0.0, 0.0), plasmid_supported=False, in_provirus=True), ("chromosomal_mge", "provirus")),
])
def test_classify_plasmid_evidence(kw, expected):
    assert ac.classify(**kw) == expected


def test_main_genomad_plasmids_requires_hallmark(tmp_path):
    abr = write(tmp_path / "abr.tsv", ABR_HEAD + abr_row("S1__a", 3000, 4000, "tet(W)") + abr_row("S1__b", 3000, 4000, "tet(W)"))
    lengths = write(tmp_path / "len.tsv", "S1__a\t9000\nS1__b\t9000\n")
    scores = write(tmp_path / "agg.tsv", "seq_name\tchromosome_score\tplasmid_score\tvirus_score\n"
                   "S1__a\t0.02\t0.95\t0.03\nS1__b\t0.02\t0.95\t0.03\n")
    plas = write(tmp_path / "plas.tsv", "seq_name\tn_hallmarks\tconjugation_genes\nS1__a\t2\tNA\nS1__b\t0\tNA\n")
    out = tmp_path / "loci.tsv"
    assert ac.main(["--route", "reads", "--abricate", abr, "--lengths", lengths, "--genomad-scores", scores,
                    "--genomad-plasmids", plas, "--out", str(out)]) == 0
    rows = {r["seq"]: (r["category"], r["reason"]) for r in csv.DictReader(open(out), delimiter="\t")}
    assert rows == {"S1__a": ("plasmid", "genomad"), "S1__b": ("unresolved", "plasmid_score_only")}


def test_genes_table_is_plasmid_evidence(tmp_path):
    genes = write(tmp_path / "genes.tsv",
                  "gene\tstart\tend\tplasmid_hallmark\tannotation_conjscan\tannotation_description\n"
                  "S1__a_1\t10\t500\t1\tNA\trep protein\n"
                  "S1__b_1\t10\t500\t0\tNA\thypothetical\n"
                  "S1__c_1\t10\t500\t0\tMOBQ\trelaxase\n")
    assert ac.read_plasmid_evidence(genes) == {"S1__a", "S1__c"}
    summary = write(tmp_path / "plas.tsv", "seq_name\tn_hallmarks\tconjugation_genes\nS1__a\t2\tNA\nS1__b\t0\tNA\n")
    assert ac.read_plasmid_evidence(summary) == {"S1__a"}


def test_plasmidfinder_hit_on_the_read_is_replicon_evidence(tmp_path):
    abr = write(tmp_path / "abr.tsv", ABR_HEAD + abr_row("S1__r", 3000, 4000, "tet(L)", res="tetracycline") +
                abr_row("S1__r", 6000, 6500, "rep27_2_repA(pSGG1)", db="plasmidfinder", res=""))
    lengths = write(tmp_path / "len.tsv", "S1__r\t9000\n")
    scores = write(tmp_path / "agg.tsv", "seq_name\tchromosome_score\tplasmid_score\tvirus_score\nS1__r\t0.02\t0.95\t0.03\n")
    genes = write(tmp_path / "genes.tsv", "gene\tstart\tend\tplasmid_hallmark\tannotation_conjscan\tannotation_description\n"
                  "S1__r_1\t3000\t4000\t0\tNA\ttetracycline efflux\n")
    out = tmp_path / "loci.tsv"
    assert ac.main(["--route", "reads", "--abricate", abr, "--lengths", lengths, "--genomad-scores", scores,
                    "--genomad-genes", genes, "--genomad-plasmids", genes, "--out", str(out)]) == 0
    row = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert (row["category"], row["reason"], row["mob_detail"]) == ("plasmid", "genomad+mobsuite", "rep27_2_repA(pSGG1)")


def test_gc_delta_from_fasta(tmp_path):
    seq = "A" * 1000 + "GC" * 500 + "A" * 1000          # gene 1001-2000 is all GC
    fa = write(tmp_path / "r.fa", ">S1__r\n" + seq + "\n")
    abr = write(tmp_path / "abr.tsv", ABR_HEAD + abr_row("S1__r", 1001, 2000, "tet(W)", res="tetracycline"))
    lengths = write(tmp_path / "len.tsv", "S1__r\t3000\n")
    out = tmp_path / "loci.tsv"
    assert ac.main(["--route", "reads", "--abricate", abr, "--lengths", lengths, "--fasta", fa, "--out", str(out)]) == 0
    row = list(csv.DictReader(open(out), delimiter="\t"))[0]
    assert (row["gc_gene"], row["gc_seq"], row["gc_delta"]) == ("1.000", "0.333", "0.667")
    assert list(row) == ac.LOCUS_COLUMNS
