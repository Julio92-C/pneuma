"""extend_context.py: contig calls carried back to read hits by position."""
import csv

from conftest import load_bin

ec = load_bin("extend_context")

LOCI_HEAD = "sample\troute\tseq\tseq_len\tgene\tdrug_class\tstart\tend\tcategory\treason\n"
VARIANT_GENE = "Bifidobacterium_adolescentis_rpoB_mutants_conferring_resistance_to_rifampicin"
READ_LOCI = LOCI_HEAD + (
    "S1\treads\tS1__r1\t9000\ttet(L)\ttetracycline\t3001\t4000\tunresolved\tshort_flank\n"
    "S1\treads\tS1__r2\t9000\t%s\trifamycin\t100\t900\tchromosome\tgenomad\n"
    "S1\treads\tS1__r3\t9000\ttet(W)\ttetracycline\t100\t900\tplasmid\tgenomad\n" % VARIANT_GENE)
CONTIG_LOCI = LOCI_HEAD + (
    "S1\tcontigs\tS1__contig_1\t50000\ttet(L)\ttetracycline\t4101\t5100\tplasmid\tmobsuite\n"
    "S1\tcontigs\tS1__contig_9\t50000\ttet(W)\ttetracycline\t100\t900\tunresolved\tlow_score\n")


def write(path, text):
    path.write_text(text)
    return str(path)


def rows(text):
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def paf(q, qs, qe, strand, t, ts, te, matches):
    return "\t".join(map(str, [q, 9000, qs, qe, strand, t, 50000, ts, te, matches, qe - qs, 60])) + "\n"


def test_best_alignment_is_keyed_by_the_contig_sample(tmp_path):
    p = write(tmp_path / "a.paf", paf("r1", 0, 9000, "+", "S1__contig_1", 1000, 10000, 8500)
              + paf("r1", 0, 5000, "+", "S1__contig_2", 0, 5000, 4000))
    assert ec.read_best_alignments([p]) == {"S1__r1": ("S1__contig_1", 1000, 10000, 0, "+", 8500)}


def test_project_on_both_strands():
    assert ec.project(3001, 4000, ("c", 1000, 10000, 0, 9000, "+", 1)) == (4001, 5000)
    assert ec.project(3001, 4000, ("c", 1000, 10000, 0, 9000, "-", 1)) == (6001, 7000)
    # worked 1-based examples: read 3001-4000 (block [0,9000) -> contig [1000,10000))
    assert ec.project(1, 1, ("c", 1000, 10000, 0, 9000, "+", 1)) == (1001, 1001)
    assert ec.project(1, 1, ("c", 1000, 10000, 0, 9000, "-", 1)) == (10000, 10000)
    assert ec.project(9000, 9000, ("c", 1000, 10000, 0, 9000, "-", 1)) == (1001, 1001)


def test_read_alignments_keeps_every_alignment_with_its_query_block(tmp_path):
    p = write(tmp_path / "a.paf", paf("r1", 0, 9000, "+", "S1__contig_1", 1000, 10000, 8500)
              + paf("r1", 0, 5000, "-", "S1__contig_2", 0, 5000, 4000))
    assert ec.read_alignments([p]) == {"S1__r1": [("S1__contig_1", 1000, 10000, 0, 9000, "+", 8500, 9000, ""),
                                                  ("S1__contig_2", 0, 5000, 0, 5000, "-", 4000, 5000, "")]}


def alns(tmp_path, text):
    return ec.read_alignments([write(tmp_path / "x.paf", text)])


def test_extend_takes_a_resolved_contig_call_at_the_projected_position(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_1", 1000, 10000, 8500)
              + paf("r3", 0, 9000, "+", "S1__contig_9", 0, 9000, 8500))
    by = {r["seq"]: r for r in ec.extend(rows(READ_LOCI), rows(CONTIG_LOCI), al, "card", pad=500)}
    assert (by["S1__r1"]["contig"], by["S1__r1"]["contig_context"]) == ("S1__contig_1", "plasmid")
    assert (by["S1__r1"]["final_context"], by["S1__r1"]["context_source"]) == ("plasmid", "contig")
    assert (by["S1__r2"]["contig_context"], by["S1__r2"]["final_context"]) == ("no_contig_match", "chromosome")
    # an unresolved contig never overrides a resolved read call
    assert by["S1__r3"]["contig_context"] == "unresolved"
    assert (by["S1__r3"]["final_context"], by["S1__r3"]["context_source"]) == ("plasmid", "read")


def one_read(gene, start, end):
    return rows(LOCI_HEAD + "S1\treads\tS1__r1\t9000\t%s\tx\t%d\t%d\tunresolved\tshort_flank\n" % (gene, start, end))


def contig_hits(*hits):
    return rows(LOCI_HEAD + "".join("S1\tcontigs\t%s\t50000\t%s\tx\t%d\t%d\t%s\tr\n" % h for h in hits))


def test_extend_minus_strand_alignment(tmp_path):
    # gene at read 3001-4000, read on the minus strand of contig 1000-10000 -> contig 6000-6999
    al = alns(tmp_path, paf("r1", 0, 9000, "-", "S1__contig_2", 1000, 10000, 8500))
    hits = contig_hits(("S1__contig_2", "tet(L)", 4001, 5000, "chromosome"),
                       ("S1__contig_2", "tet(L)", 6001, 7000, "plasmid"))
    (r,) = ec.extend(one_read("tet(L)", 3001, 4000), hits, al, "card")
    assert (r["contig"], r["contig_context"], r["final_context"]) == ("S1__contig_2", "plasmid", "plasmid")


def test_extend_uses_the_alignment_that_covers_the_gene_not_the_best_one(tmp_path):
    # best alignment (most matches) covers read 0-6000; the gene (7001-8000) is only in the supplementary block
    al = alns(tmp_path, paf("r1", 0, 6000, "+", "S1__contig_3", 0, 6000, 5900)
              + paf("r1", 6500, 9000, "+", "S1__contig_4", 100, 2600, 2400))
    hits = contig_hits(("S1__contig_3", "tet(L)", 7001, 8000, "chromosome"),
                       ("S1__contig_4", "tet(L)", 601, 1600, "plasmid"))
    (r,) = ec.extend(one_read("tet(L)", 7001, 8000), hits, al, "card")
    assert (r["contig"], r["contig_context"], r["context_source"]) == ("S1__contig_4", "plasmid", "contig")


def test_extend_no_alignment_covering_the_gene_is_no_contig_match(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 6000, "+", "S1__contig_3", 0, 6000, 5900))
    hits = contig_hits(("S1__contig_3", "tet(L)", 7001, 8000, "plasmid"))
    (r,) = ec.extend(one_read("tet(L)", 7001, 8000), hits, al, "card")
    assert (r["contig_context"], r["final_context"], r["context_source"]) == ("no_contig_match", "unresolved", "read")


def test_extend_tandem_genes_same_gene_wins_over_a_nearer_one(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 8800))
    hits = contig_hits(("S1__contig_5", "qacEdelta1", 1001, 1800, "plasmid"),
                       ("S1__contig_5", "sul1", 1201, 2000, "chromosomal_mge"))
    (r,) = ec.extend(one_read("sul1", 1001, 1800), hits, al, "card")
    assert (r["contig_gene"], r["contig_context"]) == ("sul1", "chromosomal_mge")


def test_extend_larger_overlap_wins_among_other_genes(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 8800))
    hits = contig_hits(("S1__contig_5", "geneA", 1601, 2400, "chromosome"),
                       ("S1__contig_5", "geneB", 1101, 1900, "plasmid"))
    (r,) = ec.extend(one_read("sul1", 1001, 1800), hits, al, "card")
    assert (r["contig_gene"], r["contig_context"]) == ("geneB", "plasmid")


def test_extend_resolved_beats_unresolved_at_equal_overlap(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_6", 0, 9000, 8800))
    hits = contig_hits(("S1__contig_6", "tet(L)", 1001, 2000, "unresolved"),
                       ("S1__contig_6", "tet(L)", 1001, 2000, "plasmid"))
    (r,) = ec.extend(one_read("tet(L)", 1001, 2000), hits, al, "card")
    assert (r["contig_context"], r["final_context"], r["context_source"]) == ("plasmid", "plasmid", "contig")


def test_card_variant_models_are_flagged_only_for_card():
    flags = {r["seq"]: r["intrinsic_variant"] for r in ec.extend(rows(READ_LOCI), [], {}, "card")}
    assert flags == {"S1__r1": "no", "S1__r2": "yes", "S1__r3": "no"}
    assert {r["intrinsic_variant"] for r in ec.extend(rows(READ_LOCI), [], {}, "vfdb")} == {"no"}


def test_cli_with_no_contigs_writes_the_header(tmp_path):
    rl, cl, out = write(tmp_path / "r.tsv", LOCI_HEAD), write(tmp_path / "c.tsv", LOCI_HEAD), tmp_path / "o.tsv"
    assert ec.main(["--db", "card", "--read-loci", rl, "--contig-loci", cl, "--paf", "--out", str(out)]) == 0
    assert out.read_text().splitlines() == ["\t".join(ec.COLUMNS)]


def test_extend_weak_alignment_over_the_gene_is_no_contig_match(tmp_path):
    # covers the gene but identity 5000 / 9000 = 0.56 < 0.8: ignored, the read call stays
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 5000))
    hits = contig_hits(("S1__contig_5", "tet(L)", 1001, 1800, "plasmid"))
    (r,) = ec.extend(one_read("tet(L)", 1001, 1800), hits, al, "card")
    assert (r["contig_context"], r["final_context"], r["context_source"]) == ("no_contig_match", "unresolved", "read")
    # the floor is a parameter
    (r,) = ec.extend(one_read("tet(L)", 1001, 1800), hits, al, "card", min_identity=0.5)
    assert (r["contig_context"], r["context_source"]) == ("plasmid", "contig")


def test_extend_weak_best_alignment_is_skipped_for_a_passing_one(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 6000)
              + paf("r1", 500, 2500, "+", "S1__contig_7", 0, 2000, 1900))
    hits = contig_hits(("S1__contig_5", "tet(L)", 1001, 1800, "chromosome"),
                       ("S1__contig_7", "tet(L)", 501, 1300, "plasmid"))
    (r,) = ec.extend(one_read("tet(L)", 1001, 1800), hits, al, "card")
    assert (r["contig"], r["contig_context"]) == ("S1__contig_7", "plasmid")


def test_extend_different_gene_near_but_not_overlapping_is_no_match(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 8800))
    hits = contig_hits(("S1__contig_5", "qacEdelta1", 1901, 2300, "plasmid"))   # 100 bp past the gene, within pad
    (r,) = ec.extend(one_read("sul1", 1001, 1800), hits, al, "card")
    assert (r["contig_context"], r["final_context"], r["context_source"]) == ("no_contig_match", "unresolved", "read")


def test_extend_same_gene_near_but_not_overlapping_matches(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 8800))
    hits = contig_hits(("S1__contig_5", "sul1", 1901, 2700, "plasmid"))
    (r,) = ec.extend(one_read("sul1", 1001, 1800), hits, al, "card")
    assert (r["contig_gene"], r["contig_context"]) == ("sul1", "plasmid")


def test_extend_different_gene_overlapping_the_projection_matches(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 8800))
    hits = contig_hits(("S1__contig_5", "qacEdelta1", 1701, 2100, "plasmid"))
    (r,) = ec.extend(one_read("sul1", 1001, 1800), hits, al, "card")
    assert (r["contig_gene"], r["contig_context"], r["context_source"]) == ("qacEdelta1", "plasmid", "contig")


def test_cli_min_identity_flag(tmp_path):
    rl = write(tmp_path / "r.tsv", LOCI_HEAD + "S1\treads\tS1__r1\t9000\ttet(L)\tx\t1001\t1800\tunresolved\tf\n")
    cl = write(tmp_path / "c.tsv", LOCI_HEAD + "S1\tcontigs\tS1__contig_5\t50000\ttet(L)\tx\t1001\t1800\tplasmid\tr\n")
    p = write(tmp_path / "a.paf", paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 6300))   # identity 0.70
    for flag, want in ((["--min-identity", "0.6"], "contig"), ([], "read")):
        out = tmp_path / "o.tsv"
        args = ["--db", "card", "--read-loci", rl, "--contig-loci", cl, "--paf", p, "--out", str(out)]
        assert ec.main(args + flag) == 0
        assert rows(out.read_text())[0]["context_source"] == want


def test_paf_cigar_tag_is_kept_other_tags_ignored(tmp_path):
    # minimap2 -c adds tags (NM, cg, ...) after column 12; only the cg CIGAR is kept
    line = paf("r1", 0, 9000, "+", "S1__contig_1", 1000, 10000, 8500).rstrip("\n") + "\ttp:A:P\tNM:i:500\tcg:Z:9000M\n"
    p = write(tmp_path / "a.paf", line)
    assert ec.read_alignments([p]) == {"S1__r1": [("S1__contig_1", 1000, 10000, 0, 9000, "+", 8500, 9000, "9000M")]}
    assert ec.read_best_alignments([p]) == {"S1__r1": ("S1__contig_1", 1000, 10000, 0, "+", 8500)}


def cg(contig, ts, te, strand, cigar, qs=0, qe=9000):
    return (contig, ts, te, qs, qe, strand, 8900, qe - qs, cigar)


def test_project_through_a_deletion_on_the_plus_strand():
    # 10 contig bases missing from the read after read base 3000: the gene image moves 10 bp right
    aln = cg("c", 1000, 10010, "+", "3000M10D6000M")
    assert ec.project(3001, 4000, aln) == (4011, 5010)
    assert ec.project(3001, 4000, aln[:8] + ("",)) == (4001, 5000)        # no CIGAR: linear


def test_project_through_an_insertion_on_the_plus_strand():
    # read bases 3001-3010 are inserted: the gene start falls in the insertion -> next contig base
    assert ec.project(3001, 4000, cg("c", 1000, 9990, "+", "3000M10I5990M")) == (4001, 4990)


def test_project_minus_strand_all_match_equals_linear():
    aln = cg("c", 1000, 10000, "-", "9000M")
    assert ec.project(3001, 4000, aln) == (6001, 7000)
    assert ec.project(3001, 4000, aln[:8] + ("",)) == (6001, 7000)


def test_project_minus_strand_deletion_beyond_the_gene():
    # the CIGAR runs along the contig: the first 6000 bases (read 3001-9000 reversed) carry the gene;
    # a 10 bp deletion after them shifts only the linear estimate
    aln = cg("c", 1000, 10010, "-", "6000M10D3000M")
    assert ec.project(3001, 4000, aln) == (6001, 7000)
    assert ec.project(3001, 4000, aln[:8] + ("",)) == (6011, 7010)


def test_project_outside_the_block_extends_linearly():
    aln = cg("c", 1000, 10000, "+", "9000M", qs=100, qe=9100)
    assert ec.project(50, 60, aln) == (950, 960)


def test_a_cigar_that_does_not_span_the_block_falls_back_to_linear():
    aln = cg("c", 1000, 10010, "+", "100M")                               # 100 query bases for a 9000 bp block
    assert ec.project(3001, 4000, aln) == (4001, 5000)


def test_extend_uses_the_cigar_projection(tmp_path):
    # a different gene overlaps the CIGAR image (4011-5010) but not the linear one (4001-5000)
    line = paf("r1", 0, 9000, "+", "S1__contig_8", 1000, 10010, 8900).rstrip("\n") + "\tcg:Z:3000M10D6000M\n"
    hits = contig_hits(("S1__contig_8", "qacEdelta1", 5005, 5600, "plasmid"))
    (r,) = ec.extend(one_read("tet(L)", 3001, 4000), hits, alns(tmp_path, line), "card")
    assert (r["contig_gene"], r["context_source"]) == ("qacEdelta1", "contig")
    plain = paf("r1", 0, 9000, "+", "S1__contig_8", 1000, 10010, 8900)
    (r,) = ec.extend(one_read("tet(L)", 3001, 4000), hits, alns(tmp_path, plain), "card")
    assert r["contig_context"] == "no_contig_match"


def test_context_conflict_when_a_contig_call_replaces_a_different_resolved_read_call(tmp_path):
    al = alns(tmp_path, paf("r1", 0, 9000, "+", "S1__contig_5", 0, 9000, 8800))
    read = rows(LOCI_HEAD + "S1\treads\tS1__r1\t9000\ttet(L)\tx\t1001\t1800\tchromosome\tgenomad\n")
    (r,) = ec.extend(read, contig_hits(("S1__contig_5", "tet(L)", 1001, 1800, "plasmid")), al, "card")
    assert (r["final_context"], r["context_conflict"]) == ("plasmid", "yes")
    (r,) = ec.extend(one_read("tet(L)", 1001, 1800), contig_hits(("S1__contig_5", "tet(L)", 1001, 1800, "plasmid")),
                     al, "card")
    assert (r["final_context"], r["context_conflict"]) == ("plasmid", "no")           # read was unresolved


def test_gc_delta_is_carried_from_the_read_call():
    read = rows("sample\troute\tseq\tseq_len\tgene\tdrug_class\tstart\tend\tcategory\treason\tgc_delta\n"
                "S1\treads\tS1__r1\t9000\ttet(L)\tx\t1001\t1800\tunresolved\tf\t0.042\n")
    (r,) = ec.extend(read, [], {}, "card")
    assert r["gc_delta"] == "0.042"
    assert ec.extend(one_read("tet(L)", 1001, 1800), [], {}, "card")[0]["gc_delta"] == ""
