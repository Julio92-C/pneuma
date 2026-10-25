#!/usr/bin/env python3
"""Coverage gate for one reconstruction unit, from `samtools coverage` output.

  coverage_gate.py --unit <id> --coverage <id>.coverage.txt --min-depth 10 --min-breadth 0.5 \
      --out <id>.gate.tsv

Genome-wide values over all reference sequences (chromosome + plasmids):
  mean_depth = sum(meandepth_i * length_i) / sum(length_i)
  breadth    = sum(covbases_i) / sum(length_i)    (fraction of reference bases covered >= 1x)
pass = mean_depth >= min_depth and breadth >= min_breadth.
"""
import argparse
import csv
import sys

HEADER = ["#rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
OUT_COLUMNS = ["unit_id", "reads_mapped", "mean_depth", "breadth", "gate_min_depth", "gate_min_breadth", "pass"]


def gate(unit, coverage, min_depth, min_breadth):
    with open(coverage, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows or rows[0][:len(HEADER)] != HEADER:
        raise ValueError("%s is not samtools coverage output (header %s)" % (coverage, rows[0] if rows else "missing"))
    length = covered = reads = 0
    depth_sum = 0.0
    for r in rows[1:]:
        if len(r) < len(HEADER):
            continue
        n = int(r[2]) - int(r[1]) + 1
        length += n
        covered += int(r[4])
        reads += int(r[3])
        depth_sum += float(r[6]) * n
    depth = depth_sum / length if length else 0.0
    breadth = covered / length if length else 0.0
    return {"unit_id": unit, "reads_mapped": reads, "mean_depth": "%.3f" % depth, "breadth": "%.4f" % breadth,
            "gate_min_depth": min_depth, "gate_min_breadth": min_breadth,
            "pass": "true" if (depth >= min_depth and breadth >= min_breadth and reads > 0) else "false"}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--unit", required=True)
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--min-depth", type=float, required=True)
    ap.add_argument("--min-breadth", type=float, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        row = gate(args.unit, args.coverage, args.min_depth, args.min_breadth)
    except ValueError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
