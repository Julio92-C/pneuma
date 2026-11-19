#!/usr/bin/env python3
"""Shares of the final genomic context of the read hits per sample, drug class, metadata group and gene.

Input: context/hit_context.tsv (taxon_gene_context.py); each hit counts once. n_intrinsic / intrinsic_share: the hits
of the row that are CARD protein-variant / mutation models (intrinsic_variant = yes: the gene is found but its
resistance mutation is not checked), so they can be set aside.
Writes in --outdir
  arg_context_by_sample.tsv  sample, database, category, n_hits, fraction, n_intrinsic, intrinsic_share
  arg_context_by_class.tsv   CARD hits per drug class (a multi-class gene counts in each): drug_class, category, ...
  gene_context.tsv           database, gene, n_hits, n_samples, hits per category, n_intrinsic
  arg_context_by_group.tsv   with --metadata and --group-cols: group_by, group, database, category, n_samples, n_hits,
                             pooled_fraction (hits summed over the group), mean_sample_fraction, n_intrinsic
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

CATEGORIES = ["plasmid", "chromosomal_mge", "chromosome", "unresolved"]
SAMPLE_COLUMNS = ["sample", "database", "category", "n_hits", "fraction", "n_intrinsic", "intrinsic_share"]
CLASS_COLUMNS = ["drug_class", "category", "n_hits", "fraction", "n_intrinsic", "intrinsic_share"]
GENE_COLUMNS = ["database", "gene", "n_hits", "n_samples"] + CATEGORIES + ["n_intrinsic"]
GROUP_COLUMNS = ["group_by", "group", "database", "category", "n_samples", "n_hits", "pooled_fraction",
                 "mean_sample_fraction", "n_intrinsic"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def ratio(a, b):
    return "%.4f" % (a / b) if b else "NA"


def shares(rows, keys):
    """One row per key tuple and category (all four): hits, the category's share of the key's hits, intrinsic hits."""
    n, intr, tot = defaultdict(int), defaultdict(int), defaultdict(int)
    for r in rows:
        k = tuple(r[x] for x in keys)
        n[k + (r["final_context"],)] += 1
        intr[k + (r["final_context"],)] += r["intrinsic_variant"] == "yes"
        tot[k] += 1
    out = []
    for k in sorted(tot):
        for c in CATEGORIES:
            m, i = n[k + (c,)], intr[k + (c,)]
            out.append(dict(zip(keys, k), category=c, n_hits=m, fraction=ratio(m, tot[k]), n_intrinsic=i,
                            intrinsic_share=ratio(i, m)))
    return out


def by_class(rows):
    split = []
    for r in rows:
        if r["database"] != "card":
            continue
        classes = [c.strip() for c in (r.get("drug_class") or "").split(";") if c.strip()] or ["unknown"]
        split.extend(dict(r, drug_class=c) for c in classes)
    return shares(split, ["drug_class"])


def gene_context(rows):
    acc = {}
    for r in rows:
        g = acc.setdefault((r["database"], r["gene"]), dict({c: 0 for c in CATEGORIES}, n_hits=0, samples=set(),
                                                             n_intrinsic=0))
        g["n_hits"] += 1
        g["samples"].add(r["sample"])
        g[r["final_context"]] += 1
        g["n_intrinsic"] += r["intrinsic_variant"] == "yes"
    return [dict({c: g[c] for c in CATEGORIES}, database=db, gene=gene, n_hits=g["n_hits"],
                 n_samples=len(g["samples"]), n_intrinsic=g["n_intrinsic"]) for (db, gene), g in sorted(acc.items())]


def by_group(rows, meta, group_cols):
    per = defaultdict(lambda: defaultdict(int))     # (sample, db) -> category -> hits
    intr = defaultdict(int)                         # (sample, db, category) -> intrinsic hits
    for r in rows:
        per[(r["sample"], r["database"])][r["final_context"]] += 1
        intr[(r["sample"], r["database"], r["final_context"])] += r["intrinsic_variant"] == "yes"
    out = []
    for col in group_cols:
        cells = defaultdict(list)
        for sample, db in per:
            g = (meta.get(sample, {}).get(col) or "").strip()
            if g:
                cells[(g, db)].append(sample)
        for g, db in sorted(cells):
            samples = sorted(cells[(g, db)])
            totals = {s: sum(per[(s, db)].values()) for s in samples}
            for c in CATEGORIES:
                hits = sum(per[(s, db)][c] for s in samples)
                fracs = [per[(s, db)][c] / totals[s] for s in samples]
                out.append({"group_by": col, "group": g, "database": db, "category": c, "n_samples": len(samples),
                            "n_hits": hits, "pooled_fraction": ratio(hits, sum(totals.values())),
                            "mean_sample_fraction": "%.4f" % (sum(fracs) / len(fracs)),
                            "n_intrinsic": sum(intr[(s, db, c)] for s in samples)})
    return out


def write(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hit-context", required=True, help="context/hit_context.tsv")
    ap.add_argument("--metadata", default="", help="CSV with a 'sample' column (--meridian_metadata)")
    ap.add_argument("--group-cols", default="", help="comma-separated metadata columns (--meridian_group)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)
    rows = read_tsv(a.hit_context)
    os.makedirs(a.outdir, exist_ok=True)
    write(os.path.join(a.outdir, "arg_context_by_sample.tsv"), SAMPLE_COLUMNS, shares(rows, ["sample", "database"]))
    write(os.path.join(a.outdir, "arg_context_by_class.tsv"), CLASS_COLUMNS, by_class(rows))
    write(os.path.join(a.outdir, "gene_context.tsv"), GENE_COLUMNS, gene_context(rows))
    cols = [c.strip() for c in a.group_cols.split(",") if c.strip()]
    if a.metadata and cols:
        with open(a.metadata, newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            fieldnames = reader.fieldnames or []
            meta = {r["sample"]: r for r in reader} if "sample" in fieldnames else {}
        missing = [c for c in ["sample"] + cols if c not in fieldnames]
        if missing:
            print("ERROR: %s has no column(s) %s" % (a.metadata, ", ".join(missing)), file=sys.stderr)
            return 1
        write(os.path.join(a.outdir, "arg_context_by_group.tsv"), GROUP_COLUMNS, by_group(rows, meta, cols))
    print("%d hits summarised" % len(rows), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
