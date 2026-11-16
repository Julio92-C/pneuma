"""circular_candidates.py: circular contigs carrying a hit become plasmid drafts."""
import csv

from conftest import load_bin

ci = load_bin("circular_candidates")

ABR = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY"
       "\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n")


def hit(seq, gene, db):
    return "c.fa\t%s\t10\t500\t+\t%s\t1-1/1\t===\t0/0\t100\t99\t%s\tX\tp\tr\n" % (seq, gene, db)


def write(path, text):
    path.write_text(text)
    return str(path)


def test_candidates_need_circular_and_a_hit(tmp_path):
    info = write(tmp_path / "S1.assembly_info.tsv", "sample\tcontig\tlength\tcoverage\tcircular\n"
                 "S1\tS1__contig_1\t4000\t40\tyes\nS1\tS1__contig_2\t9000\t8\tno\nS1\tS1__contig_3\t3000\t9\tyes\n")
    hits = write(tmp_path / "hits.tsv", ABR + hit("S1__contig_1", "tet(L)", "card")
                 + hit("S1__contig_1", "rep27_2_repA(pSGG1)", "plasmidfinder") + hit("S1__contig_2", "ErmB", "card"))
    fasta = write(tmp_path / "all.fa", ">S1__contig_1\n" + "A" * 4000 + "\n>S1__contig_2\nC\n>S1__contig_3\nG\n")
    scores = write(tmp_path / "agg.tsv", "seq_name\tchromosome_score\tplasmid_score\tvirus_score\nS1__contig_1\t0.02\t0.95\t0.03\n")
    plas = write(tmp_path / "ps.tsv", "seq_name\tlength\tn_hallmarks\tconjugation_genes\nS1__contig_1\t4000\t1\tNA\n")
    mob = write(tmp_path / "mob.txt", "sample_id\trep_type(s)\trelaxase_type(s)\tmpf_type\tpredicted_mobility\n"
                "S1__contig_1\trep_cluster_1\tMOBV\t-\tmobilizable\n")
    paf = write(tmp_path / "S1.paf", "".join(
        "\t".join(map(str, [r, 4000, 0, 3900, "+", "S1__contig_1", 4000, 0, 3900, 3800, 3900, 60])) + "\n"
        for r in ("r1", "r2", "r3")))
    taxa = write(tmp_path / "taxa.tsv", "S1__r1\t1350\nS1__r2\t1350\nS1__r3\t0\n")
    prefix = tmp_path / "circular_plasmid_candidates"
    assert ci.main(["--assembly-info", info, "--hits", hits, "--fasta", fasta, "--genomad-scores", scores,
                    "--genomad-plasmids", plas, "--mob-typer", mob, "--paf", paf, "--taxa", taxa,
                    "--out-prefix", str(prefix)]) == 0
    rows = list(csv.DictReader(open(str(prefix) + ".tsv"), delimiter="\t"))
    assert [r["contig"] for r in rows] == ["S1__contig_1"]
    r = rows[0]
    assert (r["sample"], r["length"], r["reads"], r["plasmid_score"], r["hallmarks"]) == ("S1", "4000", "3", "0.95", "1")
    assert (r["replicons"], r["ARGs"], r["relaxase"], r["mpf"], r["mobility"]) == (
        "rep27_2_repA(pSGG1)", "tet(L)", "MOBV", "", "mobilizable")
    assert r["read_taxa"] == "1350:2"
    assert (tmp_path / "per_sample" / "S1.circles.fa").read_text().startswith(">S1__contig_1\nAAAA")


def test_no_inputs_header_only(tmp_path):
    empty = write(tmp_path / "e.tsv", "")
    prefix = tmp_path / "c"
    assert ci.main(["--assembly-info", "--hits", empty, "--fasta", empty, "--genomad-scores", empty,
                    "--genomad-plasmids", empty, "--paf", "--out-prefix", str(prefix)]) == 0
    assert open(str(prefix) + ".tsv").read() == "\t".join(ci.COLUMNS) + "\n"
    assert open(str(prefix) + ".fa").read() == ""
