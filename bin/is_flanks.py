#!/usr/bin/env python3
"""IS elements next to each ARG on a contig, and composite transposons (mobilome).

For every ABRicate hit of --db (card) on a contig: the nearest ISEScan element on each side within --window bp
(distance between the gene and the element; 0 when they overlap). A composite transposon is called when the same
IS family flanks the gene on both sides in inverted orientation (opposite transposase strands; CUPID, Microbial
Genomics 2026). Wrap-around on circular contigs is not considered.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

COLUMNS = ["sample", "contig", "gene", "start", "end", "strand", "left_is", "left_family", "left_strand",
           "left_distance", "right_is", "right_family", "right_strand", "right_distance", "composite_transposon"]


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_isescan(paths):
    out = defaultdict(list)
    for p in paths:
        for r in read_tsv(p):
            if not r.get("seqID") or r["seqID"] == "seqID":
                continue
            out[r["seqID"]].append({"family": r["family"], "start": int(r["isBegin"]), "end": int(r["isEnd"]),
                                    "strand": r.get("strand", "")})
    return out


def read_hits(path, db):
    hits = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("#") or len(f) < 12 or f[11].lower() != db:
                continue
            a, b = int(f[2]), int(f[3])
            hits.append({"contig": f[1], "start": min(a, b), "end": max(a, b), "strand": f[4], "gene": f[5]})
    return hits


def closest(pairs, window):
    """[(distance, element)] -> (distance, element) nearest within the window, or None."""
    near = [(max(0, d), e) for d, e in pairs if d <= window]
    return min(near, key=lambda x: x[0]) if near else None


def side(prefix, found):
    if not found:
        return {prefix + "_is": "", prefix + "_family": "", prefix + "_strand": "", prefix + "_distance": ""}
    d, e = found
    return {prefix + "_is": "%s:%d-%d" % (e["family"], e["start"], e["end"]), prefix + "_family": e["family"],
            prefix + "_strand": e["strand"], prefix + "_distance": d}


def flanks(hits, elements, window=500):
    rows = []
    for g in hits:
        ise = elements.get(g["contig"], [])
        left = closest([(g["start"] - e["end"], e) for e in ise if e["start"] < g["start"] and e["end"] < g["end"]], window)
        right = closest([(e["start"] - g["end"], e) for e in ise if e["start"] > g["start"] and e["end"] > g["end"]], window)
        composite = bool(left and right and left[1]["family"] == right[1]["family"]
                         and left[1]["strand"] != right[1]["strand"])
        row = {"sample": g["contig"].split("__", 1)[0], "contig": g["contig"], "gene": g["gene"], "start": g["start"],
               "end": g["end"], "strand": g["strand"]}
        row.update(side("left", left))
        row.update(side("right", right))
        row["composite_transposon"] = "yes" if composite else "no"
        rows.append(row)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--isescan", nargs="*", default=[], help="ISEScan .tsv tables")
    ap.add_argument("--hits", required=True, help="ABRicate on the contigs (tab format)")
    ap.add_argument("--db", default="card")
    ap.add_argument("--window", type=int, default=500)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    rows = flanks(read_hits(a.hits, a.db.lower()), read_isescan(a.isescan), a.window)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("%d ARG hits on contigs; %d next to an IS element; %d composite transposons" % (
        len(rows), sum(bool(r["left_is"] or r["right_is"]) for r in rows),
        sum(r["composite_transposon"] == "yes" for r in rows)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
