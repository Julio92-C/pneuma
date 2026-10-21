#!/usr/bin/env python3
"""Merge per-database row files into database_manifest.tsv.

Each input file holds tab-separated rows: name, version, path (relative to the database
output directory), sha256. The output has the header `name version path sha256`, one row per
database in a fixed order (Kraken2, taxdump, Hostile, ABRicate), with `--root` prepended to
every path. Fails (exit 1) on a row with a missing field or an empty sha256/version, and when
there are no rows at all.
"""
import argparse
import sys

ORDER = ["kraken2_pluspf", "ncbi_taxdump", "hostile_index",
         "abricate_card", "abricate_vfdb", "abricate_plasmidfinder"]


def rank(name):
    return (ORDER.index(name), "") if name in ORDER else (len(ORDER), name)


def build(row_files, root):
    rows = []
    problems = []
    for fname in row_files:
        with open(fname) as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                cols = line.split("\t")
                if len(cols) != 4 or any(not c.strip() for c in cols):
                    problems.append(
                        f"{fname}:{lineno}: incomplete row (need name, version, path, sha256): {line!r}")
                    continue
                name, version, path, sha = cols
                rows.append((name, version, f"{root.rstrip('/')}/{path}", sha))
    if not rows and not problems:
        problems.append("no database rows to write")
    if problems:
        raise ValueError("\n".join(problems))
    rows.sort(key=lambda r: rank(r[0]))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="database output directory prepended to every path")
    ap.add_argument("--out", default="database_manifest.tsv")
    ap.add_argument("rows", nargs="+", help="row files")
    args = ap.parse_args(argv)
    try:
        rows = build(args.rows, args.root)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    with open(args.out, "w") as out:
        out.write("name\tversion\tpath\tsha256\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
