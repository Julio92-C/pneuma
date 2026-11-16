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
    assert ec.project(3001, 4000, ("c", 1000, 10000, 0, "+", 1)) == (4001, 5000)
    assert ec.project(3001, 4000, ("c", 1000, 10000, 0, "-", 1)) == (6000, 6999)


def test_extend_takes_a_resolved_contig_call_at_the_projected_position():
    best = {"S1__r1": ("S1__contig_1", 1000, 10000, 0, "+", 8500),
            "S1__r3": ("S1__contig_9", 0, 9000, 0, "+", 8500)}
    by = {r["seq"]: r for r in ec.extend(rows(READ_LOCI), rows(CONTIG_LOCI), best, "card", pad=500)}
    assert (by["S1__r1"]["contig"], by["S1__r1"]["contig_context"]) == ("S1__contig_1", "plasmid")
    assert (by["S1__r1"]["final_context"], by["S1__r1"]["context_source"]) == ("plasmid", "contig")
    assert (by["S1__r2"]["contig_context"], by["S1__r2"]["final_context"]) == ("no_contig_match", "chromosome")
    # an unresolved contig never overrides a resolved read call
    assert by["S1__r3"]["contig_context"] == "unresolved"
    assert (by["S1__r3"]["final_context"], by["S1__r3"]["context_source"]) == ("plasmid", "read")


def test_card_variant_models_are_flagged_only_for_card():
    flags = {r["seq"]: r["intrinsic_variant"] for r in ec.extend(rows(READ_LOCI), [], {}, "card")}
    assert flags == {"S1__r1": "no", "S1__r2": "yes", "S1__r3": "no"}
    assert {r["intrinsic_variant"] for r in ec.extend(rows(READ_LOCI), [], {}, "vfdb")} == {"no"}


def test_cli_with_no_contigs_writes_the_header(tmp_path):
    rl, cl, out = write(tmp_path / "r.tsv", LOCI_HEAD), write(tmp_path / "c.tsv", LOCI_HEAD), tmp_path / "o.tsv"
    assert ec.main(["--db", "card", "--read-loci", rl, "--contig-loci", cl, "--paf", "--out", str(out)]) == 0
    assert out.read_text().splitlines() == ["\t".join(ec.COLUMNS)]
