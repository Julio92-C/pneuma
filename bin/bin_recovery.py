#!/usr/bin/env python3
"""Community-assembly bins: quality, read-vote taxonomy, panel status, and depth/breadth in every sample of the patient.

Taxonomy is a majority vote of the Kraken2 species of the bin sample's reads that map (primary, MAPQ >= 1) to the
bin's contigs. panel_decision is the bin sample's decision in pathogen_candidates.tsv for that species, or
not_in_panel. Persistence rows come from mapping each sample's reads to all bins of its patient.
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from arg_context import read_lengths  # noqa: E402
from pathogen_report import mimag_grade  # noqa: E402

BIN_COLUMNS = ["sample", "patient", "bin", "n_contigs", "size_bp", "completeness", "contamination", "mimag_grade",
               "species", "species_vote", "genus", "n_reads", "panel_decision", "own_depth", "own_breadth"]
PERSIST_COLUMNS = ["patient", "bin", "bin_sample", "read_sample", "mean_depth", "breadth", "present"]
TAXID_RE = re.compile(r"taxid (\d+)\)")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def strip(path, suffix):
    name = os.path.basename(path)
    if not name.endswith(suffix):
        raise ValueError("%s does not end with %s" % (name, suffix))
    return name[: -len(suffix)]


def read_taxonomy(taxdump):
    parent, rank, name = {}, {}, {}
    with open(os.path.join(taxdump, "nodes.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t")
            parent[f[0]], rank[f[0]] = f[1], f[2].rstrip("\t|\n").strip()
    with open(os.path.join(taxdump, "names.dmp")) as fh:
        for line in fh:
            f = line.rstrip("\t|\n").split("\t|\t")
            if len(f) >= 4 and f[3] == "scientific name":
                name[f[0]] = f[1]
    return parent, rank, name


def ancestor_at(taxid, level, parent, rank):
    for _ in range(100):
        if taxid not in parent:
            return None
        if rank.get(taxid) == level:
            return taxid
        if parent[taxid] == taxid:
            return None
        taxid = parent[taxid]
    return None


def read_kraken(path):
    """Kraken2 per-read output -> {read_id: taxid} for classified reads (plain or --use-names taxids)."""
    out = {}
    for line in open(path):
        f = line.rstrip("\r\n").split("\t")
        if len(f) < 3 or f[0] != "C":
            continue
        m = TAXID_RE.search(f[2])
        out[f[1]] = m.group(1) if m else f[2].strip()
    return out


def read_pairs(path):
    return [tuple(f[:2]) for f in (l.rstrip("\r\n").split("\t") for l in open(path)) if len(f) >= 2]


def vote(read_ids, taxa, level, parent, rank):
    c = Counter()
    for rid in read_ids:
        t = taxa.get(rid)
        a = ancestor_at(t, level, parent, rank) if t else None
        if a:
            c[a] += 1
    if not c:
        return None, 0.0
    top, n = c.most_common(1)[0]
    return top, n / sum(c.values())


def recover(contig2bin, lengths, samples, read2contig, kraken, coverage, checkm2, candidates, taxonomy,
            min_depth=1.0, min_breadth=0.5):
    parent, rank, name = taxonomy
    bins, bin_sample = defaultdict(list), {}
    for r in contig2bin:
        bins[r["bin"]].append(r["contig"])
        bin_sample[r["bin"]] = r["sample"]
    contig_bin = {c: b for b, cs in bins.items() for c in cs}
    reads_of = defaultdict(list)                         # bin -> reads of its own sample on its contigs
    for sample, pairs in read2contig.items():
        for rid, contig in pairs:
            b = contig_bin.get(contig)
            if b and bin_sample[b] == sample:
                reads_of[b].append(rid)
    cov = defaultdict(dict)                              # bin -> read_sample -> [length, covbases, depth*length]
    for read_sample, rows in coverage.items():
        for r in rows:
            b = contig_bin.get(r["#rname"])
            if not b:
                continue
            n = int(r["endpos"]) - int(r["startpos"]) + 1
            acc = cov[b].setdefault(read_sample, [0, 0, 0.0])
            acc[0] += n
            acc[1] += int(r["covbases"])
            acc[2] += float(r["meandepth"]) * n
    quality = {r["Name"]: r for r in checkm2}
    decision = {(r["sample"], r["species"]): r["decision"] for r in candidates}
    patient = {r["sample"]: r["patient"] for r in samples}
    bin_rows, persist_rows = [], []
    for b in sorted(bins):
        s = bin_sample[b]
        sp, frac = vote(reads_of[b], kraken.get(s, {}), "species", parent, rank)
        ge, _ = vote(reads_of[b], kraken.get(s, {}), "genus", parent, rank)
        sp_name = name.get(sp, "") if sp else ""
        q = quality.get(b, {})
        own = cov[b].get(s)
        bin_rows.append({
            "sample": s, "patient": patient.get(s, ""), "bin": b, "n_contigs": len(bins[b]),
            "size_bp": sum(lengths.get(c, 0) for c in bins[b]),
            "completeness": q.get("Completeness", ""), "contamination": q.get("Contamination", ""),
            "mimag_grade": mimag_grade(q.get("Completeness"), q.get("Contamination")),
            "species": sp_name, "species_vote": "%.3f" % frac if sp else "",
            "genus": name.get(ge, "") if ge else "", "n_reads": len(reads_of[b]),
            "panel_decision": decision.get((s, sp_name), "not_in_panel") if sp_name else "",
            "own_depth": "%.2f" % (own[2] / own[0]) if own else "",
            "own_breadth": "%.3f" % (own[1] / own[0]) if own else ""})
        for rs in sorted(cov[b]):
            n, covb, dl = cov[b][rs]
            depth, breadth = dl / n, covb / n
            persist_rows.append({"patient": patient.get(s, ""), "bin": b, "bin_sample": s, "read_sample": rs,
                                 "mean_depth": "%.2f" % depth, "breadth": "%.3f" % breadth,
                                 "present": "yes" if depth >= min_depth and breadth >= min_breadth else "no"})
    return bin_rows, persist_rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--contig2bin", required=True, help="TSV: sample, bin, contig")
    ap.add_argument("--lengths", required=True, help="seqkit fx2tab -n -i -l of the contigs")
    ap.add_argument("--samples", required=True, help="pathogen_samples.tsv (sample, patient)")
    ap.add_argument("--read2contig", nargs="+", required=True, help="<sample>.read2contig.tsv")
    ap.add_argument("--kraken", nargs="+", required=True, help="<sample>.kraken2.classifiedreads.txt")
    ap.add_argument("--coverage", nargs="+", required=True, help="<read_sample>.bins_coverage.tsv")
    ap.add_argument("--checkm2", nargs="*", default=[], help="CheckM2 quality_report.tsv")
    ap.add_argument("--candidates", required=True, help="pathogen_candidates.tsv")
    ap.add_argument("--taxdump", required=True, help="directory with nodes.dmp and names.dmp")
    ap.add_argument("--min-depth", type=float, default=1.0)
    ap.add_argument("--min-breadth", type=float, default=0.5)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    bins, persist = recover(
        read_tsv(a.contig2bin), read_lengths(a.lengths), read_tsv(a.samples),
        {strip(p, ".read2contig.tsv"): read_pairs(p) for p in a.read2contig},
        {strip(p, ".kraken2.classifiedreads.txt"): read_kraken(p) for p in a.kraken},
        {strip(p, ".bins_coverage.tsv"): read_tsv(p) for p in a.coverage},
        [r for p in a.checkm2 for r in read_tsv(p)], read_tsv(a.candidates), read_taxonomy(a.taxdump),
        a.min_depth, a.min_breadth)
    write(os.path.join(a.outdir, "recovered_bins.tsv"), BIN_COLUMNS, bins)
    write(os.path.join(a.outdir, "bin_persistence.tsv"), PERSIST_COLUMNS, persist)
    print("%d bins, %d bin x sample rows" % (len(bins), len(persist)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
