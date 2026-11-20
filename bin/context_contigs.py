#!/usr/bin/env python3
"""Rename one sample's Flye contigs "<sample>__contig_N" and tabulate them (genomic context, layer 2).

Writes <prefix>.contigs.fa and <prefix>.assembly_info.tsv (sample, contig, length, coverage, circular yes/no, repeat yes/no: Flye's repeat flag) from
Flye's assembly.fasta(.gz) and assembly_info.txt. The sample prefix lets contigs of all samples be pooled.
"""
import argparse
import csv
import gzip
import sys

INFO_COLUMNS = ["sample", "contig", "length", "coverage", "circular", "repeat"]


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    """-> [(name, sequence)] in file order."""
    records, name, seq = [], None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def read_info(path):
    with open(path, newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r.get("#seq_name")]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fasta", required=True, help="Flye assembly.fasta(.gz)")
    ap.add_argument("--info", required=True, help="Flye assembly_info.txt")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args(argv)
    prefix = a.sample + "__"
    records = read_fasta(a.fasta)
    with open(a.out_prefix + ".contigs.fa", "w", newline="\n") as fh:
        for name, seq in records:
            fh.write(">%s%s\n%s\n" % (prefix, name, seq))
    with open(a.out_prefix + ".assembly_info.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(INFO_COLUMNS)
        for r in read_info(a.info):
            w.writerow([a.sample, prefix + r["#seq_name"], r["length"], r["cov."], "yes" if r["circ."] == "Y" else "no",
                        "yes" if r.get("repeat") == "Y" else "no"])
    print("%s: %d contigs" % (a.sample, len(records)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
