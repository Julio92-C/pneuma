#!/usr/bin/env python3
"""Carry the genomic context of locally assembled contigs back to each read hit (genomic context, layer 2).

A read hit takes the call of a contig hit of the same database at the gene's projected position: of the read's
alignments to its own sample's contigs with identity (matches / block length) >= --min-identity (0.8), the one whose
query block overlaps the gene most (ties: more matches) carries the gene's read coordinates onto the contig
(strand-aware); no such alignment over the gene means no contig match. A contig hit of the same gene matches within
--pad bp of the projected span; a hit of a different gene only when it overlaps the span. Among matching hits, the
same gene wins, then the largest overlap (or the nearest), then a resolved call over an unresolved one. The
contig call becomes the final call when it is resolved; otherwise (no contig match, or an unresolved contig) the read
call is kept. CARD protein-variant / mutation models are flagged intrinsic_variant = yes: ABRicate finds the gene but
cannot check the resistance mutation.
"""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict

RESOLVED = ("plasmid", "chromosomal_mge", "chromosome")
CATEGORIES = list(RESOLVED) + ["unresolved"]
COLUMNS = ["sample", "seq", "database", "gene", "drug_class", "read_len", "start", "end", "read_context",
           "read_reason", "contig", "contig_len", "contig_gene", "contig_context", "contig_reason",
           "final_context", "context_source", "intrinsic_variant"]
VARIANT = re.compile(r"mutant|mutation|conferring_resistance", re.IGNORECASE)


def sample_of(seq):
    return seq.split("__", 1)[0]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_best_alignments(paths):
    """PAFs of recruited reads vs their sample's contigs -> {"<sample>__<read>": (contig, ts, te, qs, strand, matches)},
    the alignment with the most matches per read (read counts and abundance in circular_candidates.py)."""
    best = {}
    for path in paths:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                q, qs, strand, t, ts, te, matches = f[0], int(f[2]), f[4], f[5], int(f[7]), int(f[8]), int(f[9])
                key = sample_of(t) + "__" + q
                if key not in best or matches > best[key][5]:
                    best[key] = (t, ts, te, qs, strand, matches)
    return best


def read_alignments(paths):
    """PAFs of recruited reads vs their sample's contigs -> {"<sample>__<read>": [(contig, ts, te, qs, qe, strand,
    matches, block), ...]}, every alignment of the read (primary, secondary and supplementary) in file order."""
    alns = defaultdict(list)
    for path in paths:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                alns[sample_of(f[5]) + "__" + f[0]].append(
                    (f[5], int(f[7]), int(f[8]), int(f[2]), int(f[3]), f[4], int(f[9]), int(f[10])))
    return dict(alns)


def overlap(a0, a1, b0, b1):
    """Overlap (bp) of two inclusive spans; <= 0 is minus the gap between them."""
    return min(a1, b1) - max(a0, b0) + 1


def identity(aln):
    """Residue matches / alignment block length (PAF columns 10 / 11)."""
    return aln[6] / aln[7] if aln[7] else 0.0


def gene_alignment(start, end, alns, min_identity=0.8):
    """Of the alignments with identity >= min_identity, the one whose query block [qs, qe) overlaps the gene
    (1-based start..end) most, ties on matches."""
    best, key = None, None
    for aln in alns:
        if identity(aln) < min_identity:
            continue
        ov = min(aln[4], end) - max(aln[3], start - 1)
        if ov > 0 and (key is None or (ov, aln[6]) > key):
            best, key = aln, (ov, aln[6])
    return best


def project(start, end, aln):
    """Gene span on the read -> span on the contig through the alignment (start of the aligned block as anchor)."""
    ts, te, qs, strand = aln[1], aln[2], aln[3], aln[5]
    a, b = start - qs, end - qs
    return (ts + a, ts + b) if strand == "+" else (te - b, te - a)


def contig_match(gene, lo, hi, hits, pad):
    """Contig hit for a gene projected to lo..hi: the same gene within pad bp, a different gene only when it overlaps
    lo..hi; then same gene first, largest overlap / nearest, resolved first."""
    near = [c for c in hits if overlap(int(c["start"]), int(c["end"]), lo - (0 if c["gene"] != gene else pad),
                                       hi + (0 if c["gene"] != gene else pad)) > 0]
    if not near:
        return None
    return min(near, key=lambda c: (c["gene"] != gene, -overlap(int(c["start"]), int(c["end"]), lo, hi),
                                    c["category"] not in RESOLVED))


def extend(read_rows, contig_rows, alignments, db, pad=500, min_identity=0.8):
    by_contig = defaultdict(list)
    for c in contig_rows:
        by_contig[c["seq"]].append(c)
    out = []
    for r in read_rows:
        hit = {"contig": "", "contig_len": "", "contig_gene": "", "contig_context": "no_contig_match",
               "contig_reason": ""}
        aln = gene_alignment(int(r["start"]), int(r["end"]), alignments.get(r["seq"], []), min_identity)
        if aln:
            lo, hi = project(int(r["start"]), int(r["end"]), aln)
            c = None
            if overlap(lo, hi, aln[1] - pad, aln[2] + pad) > 0:
                c = contig_match(r["gene"], lo, hi, by_contig.get(aln[0], []), pad)
            if c:
                hit = {"contig": c["seq"], "contig_len": c["seq_len"], "contig_gene": c["gene"],
                       "contig_context": c["category"], "contig_reason": c.get("reason", "")}
        from_contig = hit["contig_context"] in RESOLVED
        row = {"sample": r["sample"], "seq": r["seq"], "database": db, "gene": r["gene"],
               "drug_class": r.get("drug_class", ""), "read_len": r["seq_len"], "start": r["start"], "end": r["end"],
               "read_context": r["category"], "read_reason": r.get("reason", "")}
        row.update(hit)
        row["final_context"] = hit["contig_context"] if from_contig else r["category"]
        row["context_source"] = "contig" if from_contig else "read"
        row["intrinsic_variant"] = "yes" if db == "card" and VARIANT.search(r["gene"]) else "no"
        out.append(row)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--db", required=True, help="ABRicate database of these hits (card, vfdb, plasmidfinder...)")
    ap.add_argument("--read-loci", required=True, help="loci_reads_<db>.tsv (arg_context.py --route reads)")
    ap.add_argument("--contig-loci", required=True, help="loci_contigs_<db>.tsv (arg_context.py --route contigs)")
    ap.add_argument("--paf", nargs="*", default=[], help="recruited reads -> their sample's contigs (minimap2 PAF)")
    ap.add_argument("--pad", type=int, default=500, help="bp tolerance around the projected gene span (same gene)")
    ap.add_argument("--min-identity", type=float, default=0.8,
                    help="minimum alignment identity (matches / block length) of a read alignment used")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    rows = extend(read_tsv(a.read_loci), read_tsv(a.contig_loci), read_alignments(a.paf), a.db.lower(), a.pad,
                  a.min_identity)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = len(rows) or 1
    on_contig = sum(r["contig_context"] != "no_contig_match" for r in rows)
    final = Counter(r["final_context"] for r in rows)
    print("%s: %d hits, %.1f%% on a contig; final %s" % (a.db, len(rows), 100.0 * on_contig / n, ", ".join(
        "%s=%d" % (c, final[c]) for c in CATEGORIES)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
