#!/usr/bin/env python3
"""Convert Kraken2 --classified-out FASTQ into taxid-tagged FASTA for ABRicate.

Input header : @<read_id> [...] kraken:taxid|<taxid> [...]
Output header: ><read_id>_kraken:taxid|<taxid>

This is the read-naming contract MERIDIAN uses to link every ABRicate hit
to the taxon of the read that carries it.
"""
import argparse
import gzip
import re
import sys

TAXID_RE = re.compile(r"kraken:taxid\|(\d+)")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def tag_records(lines):
    it = iter(lines)
    for header in it:
        header = header.rstrip("\r\n")
        if not header:
            continue
        try:
            seq = next(it).rstrip("\r\n")
            next(it)
            next(it)
        except StopIteration:
            raise ValueError(f"truncated FASTQ record after header: {header[:60]}")
        if not header.startswith("@"):
            raise ValueError(f"not a FASTQ header: {header[:60]}")
        match = TAXID_RE.search(header)
        if match is None:
            raise ValueError(f"no kraken:taxid tag in header: {header[:60]}")
        read_id = header[1:].split()[0]
        yield f">{read_id}_kraken:taxid|{match.group(1)}\n{seq}\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq", help="Kraken2 classified reads (.fastq or .fastq.gz)")
    parser.add_argument("output", help="output FASTA")
    args = parser.parse_args(argv)
    n = 0
    with open_text(args.fastq) as fh, open(args.output, "w") as out:
        for record in tag_records(fh):
            out.write(record)
            n += 1
    print(f"tagged {n} reads", file=sys.stderr)


if __name__ == "__main__":
    main()
