#!/usr/bin/env python3
"""Per-sample inputs for the genomic-context analysis: the reads carrying an ABRicate hit.

Reads the Kraken2-tagged FASTA (TAG_READS) and the ABRicate CSV (ABRICATE_READS) of one sample and writes
  <prefix>.hit_reads.fa   reads with >= 1 hit, ids "<sample>__<read id>" (the sample prefix lets samples be pooled;
                          the Kraken2 tag "_kraken:taxid|<taxid>" is removed: geNomad uses "|" in provirus names)
  <prefix>.hits.tsv       the hits in ABRicate tab format (one header), SEQUENCE renamed the same way
  <prefix>.lengths.tsv    id <TAB> length
  <prefix>.taxa.tsv       id <TAB> Kraken2 taxid
"""
import argparse
import csv
import sys

TAG = "_kraken:taxid|"
HEADER = ["#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS",
          "%COVERAGE", "%IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"]


def split_tag(name):
    """r1_kraken:taxid|562 -> ("r1", "562"); an untagged id gets an empty taxid."""
    read, _, taxid = name.partition(TAG)
    return read, taxid


def read_hits(path):
    """ABRicate --csv output (the header may repeat, once per input file) -> (header, rows)."""
    header, rows = None, []
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            if row[0].startswith("#FILE"):
                header = header or row
                continue
            rows.append(row)
    return header or HEADER, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fasta", required=True, help="Kraken2-tagged reads (TAG_READS)")
    ap.add_argument("--abricate", required=True, help="ABRicate CSV of the same reads (ABRICATE_READS)")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args(argv)
    header, rows = read_hits(a.abricate)
    hit_ids = {r[1] for r in rows}
    prefix = a.sample + "__"
    with open(a.out_prefix + ".hits.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for r in rows:
            w.writerow(r[:1] + [prefix + split_tag(r[1])[0]] + r[2:])
    n = 0
    with open(a.fasta) as fi, open(a.out_prefix + ".hit_reads.fa", "w", newline="\n") as fa, \
            open(a.out_prefix + ".lengths.tsv", "w", newline="\n") as ln, \
            open(a.out_prefix + ".taxa.tsv", "w", newline="\n") as tx:
        name, seq = None, []

        def flush():
            nonlocal n
            if name in hit_ids:
                s = "".join(seq)
                read, taxid = split_tag(name)
                fa.write(">%s%s\n%s\n" % (prefix, read, s))
                ln.write("%s%s\t%d\n" % (prefix, read, len(s)))
                tx.write("%s%s\t%s\n" % (prefix, read, taxid))
                n += 1

        for line in fi:
            if line.startswith(">"):
                flush()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        flush()
    print("%s: %d hits on %d reads" % (a.sample, len(rows), n), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
