#!/usr/bin/env python3
"""Per-sample, per-class and per-group shares of ARG loci by genomic context (both routes), and route concordance.

Weight of a locus: 1 per ARG hit on a read (route reads); mean read depth over the gene (route contigs, from
samtools bedcov), so both routes count reads that carry the gene.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

CATEGORIES = ["plasmid", "chromosomal_mge", "chromosome", "unresolved"]
SAMPLE_COLUMNS = ["sample", "route", "category", "n_loci", "weight", "fraction"]
CLASS_COLUMNS = ["route", "drug_class", "category", "n_loci", "weight", "fraction"]
GROUP_COLUMNS = ["group_by", "group", "route", "category", "n_samples", "pooled_fraction", "mean_sample_fraction"]
CONCORDANCE_COLUMNS = ["sample", "gene", "reads_category", "reads_weight", "contigs_category", "contigs_weight",
                       "agree"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_bedcov(paths):
    """samtools bedcov on a seq/start0/end/gene BED -> {(seq, start1, end): mean depth}."""
    depth = {}
    for p in paths:
        for line in open(p):
            f = line.rstrip("\r\n").split("\t")
            if len(f) < 5:
                continue
            start0, end = int(f[1]), int(f[2])
            depth[(f[0], start0 + 1, end)] = float(f[-1]) / max(end - start0, 1)
    return depth


def weigh(rows, depth):
    for r in rows:
        r["weight"] = (depth.get((r["seq"], int(r["start"]), int(r["end"])), 0.0)
                       if r["route"] == "contigs" else 1.0)
    return rows


def shares(rows, keys):
    """One row per key tuple and category: n_loci, weight, and the category's share of the key's weight."""
    acc = defaultdict(lambda: [0, 0.0])
    tot = defaultdict(float)
    for r in rows:
        k = tuple(r[x] for x in keys)
        acc[k + (r["category"],)][0] += 1
        acc[k + (r["category"],)][1] += r["weight"]
        tot[k] += r["weight"]
    out = []
    for k in sorted(tot):
        for c in CATEGORIES:
            n, w = acc.get(k + (c,), [0, 0.0])
            out.append(dict(zip(keys, k), category=c, n_loci=n, weight="%.2f" % w,
                            fraction="%.4f" % (w / tot[k]) if tot[k] else "NA"))
    return out


def by_class(rows):
    split = []
    for r in rows:
        classes = [c.strip() for c in (r["drug_class"] or "").split(";") if c.strip()] or ["unknown"]
        split.extend(dict(r, drug_class=c) for c in classes)
    return shares(split, ["route", "drug_class"])


def by_group(sample_rows, meta, group_cols):
    """Pooled (weights summed) and mean per-sample share of each category within metadata groups."""
    out = []
    for col in group_cols:
        cells = defaultdict(lambda: defaultdict(list))   # (group, route) -> category -> [(weight, fraction)]
        for r in sample_rows:
            g = meta.get(r["sample"], {}).get(col, "")
            if g == "" or r["fraction"] == "NA":
                continue
            cells[(g, r["route"])][r["category"]].append((float(r["weight"]), float(r["fraction"])))
        for g, route in sorted(cells):
            cats = cells[(g, route)]
            total = sum(w for c in CATEGORIES for w, _ in cats.get(c, []))
            n = len(cats.get(CATEGORIES[0], []))
            for c in CATEGORIES:
                vals = cats.get(c, [])
                out.append({"group_by": col, "group": g, "route": route, "category": c, "n_samples": n,
                            "pooled_fraction": "%.4f" % (sum(w for w, _ in vals) / total) if total else "NA",
                            "mean_sample_fraction": "%.4f" % (sum(f for _, f in vals) / len(vals)) if vals else "NA"})
    return out


def concordance(rows):
    """Per sample and gene: the dominant category (by weight) in each route, and whether they agree."""
    best = defaultdict(lambda: defaultdict(float))   # (sample, gene, route) -> category -> weight
    for r in rows:
        best[(r["sample"], r["gene"], r["route"])][r["category"]] += r["weight"]

    def top(k):
        cats = best.get(k)
        if not cats:
            return "", 0.0
        return max(CATEGORIES, key=lambda x: cats.get(x, 0.0)), sum(cats.values())

    out = []
    for s, g in sorted({(k[0], k[1]) for k in best}):
        rc, rw = top((s, g, "reads"))
        cc, cw = top((s, g, "contigs"))
        out.append({"sample": s, "gene": g, "reads_category": rc, "reads_weight": "%.2f" % rw,
                    "contigs_category": cc, "contigs_weight": "%.2f" % cw,
                    "agree": "NA" if not rc or not cc else ("yes" if rc == cc else "no")})
    return out


def write(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--loci", nargs="+", required=True, help="arg_context.py outputs (either route)")
    ap.add_argument("--bedcov", nargs="*", default=[], help="samtools bedcov outputs for the contig loci")
    ap.add_argument("--metadata", default="", help="CSV with a 'sample' column")
    ap.add_argument("--group-cols", default="", help="comma-separated metadata columns")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    rows = weigh([r for p in a.loci for r in read_tsv(p)], read_bedcov(a.bedcov))
    per_sample = shares(rows, ["sample", "route"])
    write(os.path.join(a.outdir, "arg_context_by_sample.tsv"), SAMPLE_COLUMNS, per_sample)
    write(os.path.join(a.outdir, "arg_context_by_class.tsv"), CLASS_COLUMNS, by_class(rows))
    write(os.path.join(a.outdir, "arg_context_concordance.tsv"), CONCORDANCE_COLUMNS, concordance(rows))
    if a.metadata and a.group_cols:
        with open(a.metadata, newline="", encoding="utf-8-sig") as fh:
            meta = {r["sample"]: r for r in csv.DictReader(fh)}
        write(os.path.join(a.outdir, "arg_context_by_group.tsv"), GROUP_COLUMNS,
              by_group(per_sample, meta, [c for c in a.group_cols.split(",") if c]))
    print("%d loci summarised" % len(rows), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
