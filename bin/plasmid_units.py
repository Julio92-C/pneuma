#!/usr/bin/env python3
"""Group circular plasmid drafts across samples into plasmid units (mobilome).

Two drafts are the same unit when their MinHash identity (canonical k-mers, bottom-s sketches, Mash distance
1 - d with d = -ln(2j / (1 + j)) / k) is >= --min-identity and the shorter is >= --min-length-ratio of the longer;
units are single-linkage clusters. k-mer sets do not depend on where Flye opened the circle or on the strand.
Abundance: reads whose best alignment is on the draft, per million reads of the sample (SELECT_READS round 1).
Writes <prefix>.tsv (one row per unit, PU001... by prevalence) and <prefix>_members.tsv (one row per draft).
"""
import argparse
import csv
import hashlib
import heapq
import math
import os
import re
import statistics
import sys
from collections import Counter, defaultdict

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
UNION_COLUMNS = ["replicons", "mob_rep", "relaxase", "mpf", "ARGs", "VFs"]
UNIT_COLUMNS = ["unit", "n_samples", "samples", "n_contigs", "representative", "length"] + UNION_COLUMNS[:4] + [
    "mobility", "ARGs", "VFs", "IS_families", "read_taxa", "mean_rpm"]
MEMBER_COLUMNS = ["unit", "contig", "sample", "length", "coverage", "reads", "rpm"]


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def kmer_hash(kmer):
    return int.from_bytes(hashlib.blake2b(kmer.encode(), digest_size=8).digest(), "big")


def sketch(seq, k=21, size=1000):
    seq = seq.upper()
    rc = seq.translate(COMP)[::-1]
    n = len(seq)
    hashes = {kmer_hash(min(seq[i:i + k], rc[n - i - k:n - i])) for i in range(n - k + 1)}
    return set(heapq.nsmallest(size, hashes))


def identity(a, b, k=21, size=1000):
    union = heapq.nsmallest(size, a | b)
    if not union:
        return 0.0
    j = sum(1 for h in union if h in a and h in b) / len(union)
    if j == 0:
        return 0.0
    return max(0.0, 1 + math.log(2 * j / (1 + j)) / k)


def clusters(records, min_identity=0.95, min_length_ratio=0.8, k=21, size=1000):
    """[(name, sequence)] -> [[names]] single-linkage clusters."""
    lengths = [len(s) for _, s in records]
    sketches = [sketch(s, k, size) for _, s in records]
    parent = list(range(len(records)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            lo, hi = sorted((lengths[i], lengths[j]))
            if hi == 0 or lo / hi < min_length_ratio:
                continue
            if identity(sketches[i], sketches[j], k, size) >= min_identity:
                parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i, (name, _s) in enumerate(records):
        groups[find(i)].append(name)
    return list(groups.values())


def split(value):
    return {x for x in re.split(r"[;,]", value or "") if x and x != "-"}


def read_totals(paths):
    totals = {}
    for p in paths:
        for r in read_tsv(p):
            if r.get("round") == "1" and (r.get("total_reads") or "").isdigit():
                totals[r["sample"]] = int(r["total_reads"])
    return totals


def read_is_families(paths):
    fam = defaultdict(set)
    for p in paths:
        for r in read_tsv(p):
            if r.get("seqID") and r.get("family") and r["seqID"] != "seqID":
                fam[r["seqID"]].add(r["family"])
    return fam


def units(cands, seqs, is_families, totals, min_identity=0.95, min_length_ratio=0.8, k=21, size=1000):
    by = {r["contig"]: r for r in cands if r["contig"] in seqs}
    groups = clusters([(c, seqs[c]) for c in sorted(by)], min_identity, min_length_ratio, k, size)
    groups.sort(key=lambda g: (-len({by[c]["sample"] for c in g}), -max(int(by[c]["length"]) for c in g), min(g)))
    rows, members = [], []
    for n, group in enumerate(groups, 1):
        uid = "PU%03d" % n
        mem = sorted((by[c] for c in group), key=lambda r: r["contig"])
        rep = max(mem, key=lambda r: (int(r["reads"] or 0), int(r["length"]), r["contig"]))
        per_sample, votes = defaultdict(float), Counter()
        for r in mem:
            total = totals.get(r["sample"])
            rpm = round(int(r["reads"] or 0) * 1e6 / total, 3) if total else ""
            members.append({"unit": uid, "contig": r["contig"], "sample": r["sample"], "length": r["length"],
                            "coverage": r["coverage"], "reads": r["reads"], "rpm": rpm})
            if rpm != "":
                per_sample[r["sample"]] += rpm
            for kv in filter(None, (r.get("read_taxa") or "").split(",")):
                taxid, count = kv.rsplit(":", 1)
                votes[taxid] += int(count)
        mobility = Counter(r["mobility"] for r in mem if r.get("mobility"))
        samples = sorted({r["sample"] for r in mem})
        row = {"unit": uid, "n_samples": len(samples), "samples": ",".join(samples), "n_contigs": len(mem),
               "representative": rep["contig"], "length": rep["length"],
               "mobility": mobility.most_common(1)[0][0] if mobility else "",
               "IS_families": ";".join(sorted(set().union(*(is_families.get(r["contig"], set()) for r in mem)))),
               "read_taxa": ",".join("%s:%d" % kv for kv in votes.most_common(3)),
               "mean_rpm": round(statistics.mean(per_sample.values()), 3) if per_sample else ""}
        for col in UNION_COLUMNS:
            row[col] = ";".join(sorted(set().union(*(split(r.get(col)) for r in mem))))
        rows.append(row)
    return rows, members


def write(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--candidates", required=True, help="circular_plasmid_candidates.tsv")
    ap.add_argument("--fasta", required=True, help="circular_plasmid_candidates.fa")
    ap.add_argument("--isescan", nargs="*", default=[])
    ap.add_argument("--recruit", nargs="*", default=[], help="SELECT_READS *.recruit.tsv (round 1 gives total reads)")
    ap.add_argument("--min-identity", type=float, default=0.95)
    ap.add_argument("--min-length-ratio", type=float, default=0.8)
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--sketch-size", type=int, default=1000)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args(argv)
    rows, members = units(read_tsv(a.candidates), read_fasta(a.fasta), read_is_families(a.isescan),
                          read_totals(a.recruit), a.min_identity, a.min_length_ratio, a.k, a.sketch_size)
    write(a.out_prefix + ".tsv", UNIT_COLUMNS, rows)
    write(a.out_prefix + "_members.tsv", MEMBER_COLUMNS, members)
    print("%d plasmid drafts -> %d plasmid units (%d in more than one sample)" % (
        len(members), len(rows), sum(r["n_samples"] > 1 for r in rows)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
