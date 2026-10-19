#!/usr/bin/env python3
"""Trim a raw Recentrifuge 2.1.2 data CSV to what MERIDIAN v1.0.0 expects.

MERIDIAN (R/02 clean_data) keeps column 1 plus every third column and strips
sample names with gsub(".*?/", "") then gsub("\..*", ""). It therefore needs
column 1 followed by exactly one cnt/una/sco triplet per sample and nothing
else. The raw file also carries derived triplets (<S>.kraken2.classifiedreads_
EXCLUSIVE_<rank>, SHARED_<rank>, SHARED_SUMMARY) and Details/Name/Rank columns,
which make the stripped names non-unique.

Output: column 1 + one triplet per sample ID (in the given order), headers
rewritten to 'kraken/<S>' (no dots). Row labels (Samples/Stats/Id) and all
taxid rows are kept unchanged; ragged rows are padded.
"""
import argparse
import csv
import re
import sys

TRIPLET = ["cnt", "una", "sco"]
_SUFFIX = re.compile(r"(\.kraken2)?\.classifiedreads(\.txt)?$")


def header_sample(cell):
    """'kraken/S1.kraken2.classifiedreads' -> 'S1'; derived columns keep their suffix."""
    name = cell.strip().replace("\\", "/").rsplit("/", 1)[-1]
    return _SUFFIX.sub("", name)


def trim(rows, samples):
    if len(rows) < 3:
        raise ValueError("Recentrifuge CSV needs at least 3 rows (Samples, Stats, Id)")
    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    header, stats = rows[0], rows[1]
    dup = sorted({s for s in samples if samples.count(s) > 1})
    if dup:
        raise ValueError(f"duplicate sample IDs: {dup}")

    keep = [0]
    new_names = {}
    for s in samples:
        hits = [j for j in range(1, width)
                if header[j].strip() and header_sample(header[j]) == s]
        good = [j for j in hits if stats[j:j + 3] == TRIPLET]
        if len(good) != 1:
            raise ValueError(
                f"sample '{s}' maps to {len(good)} cnt/una/sco triplets in the "
                f"Recentrifuge CSV (header matches: {len(hits)}); expected exactly 1")
        j = good[0]
        keep.extend([j, j + 1, j + 2])
        new_names[j] = f"kraken/{s}"

    out = []
    for i, row in enumerate(rows):
        new = [row[j] for j in keep]
        if i == 0:
            new = [new_names.get(j, "") if j else row[0] for j in keep]
        out.append(new)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rcf", required=True, help="raw pneuma.rcf.data.csv")
    p.add_argument("--samples", required=True, help="comma-separated sample IDs (controls included)")
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)
    samples = [s.strip() for s in a.samples.split(",") if s.strip()]
    if not samples:
        sys.exit("error: no sample IDs given")
    with open(a.rcf, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    try:
        out = trim(rows, samples)
    except ValueError as e:
        sys.exit(f"error: {e}")
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        csv.writer(fh, lineterminator="\n").writerows(out)


if __name__ == "__main__":
    main()
