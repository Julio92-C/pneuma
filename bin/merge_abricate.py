#!/usr/bin/env python3
"""Concatenate per-sample ABRicate --csv outputs into one summary_reportClean.csv.

Lines are copied verbatim (ABRicate does not quote fields, so re-parsing
with csv could split PRODUCT values). Repeated headers and blank lines are
dropped; one header is written first.
"""
import argparse
import sys

HEADER = ("#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,"
          "%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE")


def merge(paths, out):
    out.write(HEADER + "\n")
    rows = 0
    for path in sorted(str(p) for p in paths):
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#FILE"):
                    if line != HEADER:
                        raise ValueError(f"unexpected ABRicate header in {path}: {line}")
                    continue
                out.write(line + "\n")
                rows += 1
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args(argv)
    with open(args.out, "w") as out:
        n = merge(args.inputs, out)
    print(f"merged {n} ABRicate hits from {len(args.inputs)} files", file=sys.stderr)


if __name__ == "__main__":
    main()
