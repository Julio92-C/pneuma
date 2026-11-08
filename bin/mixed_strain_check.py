#!/usr/bin/env python3
"""Mixed-strain check of a patient-level (pooled) assembly from `samtools mpileup` text.

  mixed_strain_check.py --unit <id> --mpileup <id>.mpileup[.gz] --min-depth 8 --min-maf 0.25 \
      --min-alt 3 --max-het-per-kb 1.0 --min-evaluable 10000 --out <id>.mixed.tsv

Evaluable site: >= min-depth A/C/G/T bases (indels ignored). Het site: second most common base
with frequency >= min-maf, >= min-alt reads, seen on BOTH strands (upper case / '.' = forward,
lower case / ',' = reverse). het_per_kb = 1000 * het / evaluable.
verdict = not_evaluable (evaluable < min-evaluable) | mixed (het_per_kb > max) | single.
"""
import argparse
import csv
import gzip
import sys

OUT_COLUMNS = ["unit_id", "evaluable_sites", "het_sites", "het_per_kb", "verdict",
               "min_depth", "min_maf", "min_alt", "max_het_per_kb"]


def count_bases(ref, bases):
    """{base: [forward, reverse]} for A/C/G/T from an mpileup read-base string."""
    counts = {b: [0, 0] for b in "ACGT"}
    ref = ref.upper()
    i, n = 0, len(bases)
    while i < n:
        c = bases[i]
        if c == "^":                       # read start, followed by its mapping-quality character
            i += 2
            continue
        if c in "+-":                      # indel: [+-]<length><sequence>
            j = i + 1
            while j < n and bases[j].isdigit():
                j += 1
            i = j + (int(bases[i + 1:j]) if j > i + 1 else 0)
            continue
        if c == "." and ref in counts:
            counts[ref][0] += 1
        elif c == "," and ref in counts:
            counts[ref][1] += 1
        elif c in "ACGT":
            counts[c][0] += 1
        elif c in "acgt":
            counts[c.upper()][1] += 1
        i += 1                             # '$', '*', '#', '<', '>', N/n: ignored
    return counts


def check(lines, min_depth, min_maf, min_alt):
    evaluable = het = 0
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        counts = count_bases(f[2], f[4])
        totals = sorted(((sum(v), b) for b, v in counts.items()), reverse=True)
        depth = sum(t for t, _b in totals)
        if depth < min_depth:
            continue
        evaluable += 1
        alt_n, alt_b = totals[1]
        if alt_n >= min_alt and alt_n / depth >= min_maf and counts[alt_b][0] > 0 and counts[alt_b][1] > 0:
            het += 1
    return evaluable, het


def open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--unit", required=True)
    ap.add_argument("--mpileup", required=True)
    ap.add_argument("--min-depth", type=int, default=8)
    ap.add_argument("--min-maf", type=float, default=0.25)
    ap.add_argument("--min-alt", type=int, default=3)
    ap.add_argument("--max-het-per-kb", type=float, default=1.0)
    ap.add_argument("--min-evaluable", type=int, default=10000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    with open_text(a.mpileup) as fh:
        evaluable, het = check(fh, a.min_depth, a.min_maf, a.min_alt)
    per_kb = 1000.0 * het / evaluable if evaluable else 0.0
    if evaluable < a.min_evaluable:
        verdict = "not_evaluable"
    elif per_kb > a.max_het_per_kb:
        verdict = "mixed"
    else:
        verdict = "single"
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow({"unit_id": a.unit, "evaluable_sites": evaluable, "het_sites": het, "het_per_kb": "%.3f" % per_kb,
                    "verdict": verdict, "min_depth": a.min_depth, "min_maf": a.min_maf, "min_alt": a.min_alt,
                    "max_het_per_kb": a.max_het_per_kb})
    return 0


if __name__ == "__main__":
    sys.exit(main())
