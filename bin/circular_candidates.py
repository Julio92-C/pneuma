#!/usr/bin/env python3
"""Circular contigs carrying a genetic-element hit: draft complete plasmids (genomic context, layer 2).

Inputs
  --assembly-info     <sample>.assembly_info.tsv tables (context_contigs.py), any number
  --hits              ABRicate on all contigs (every context database, tab format, headers may repeat)
  --fasta             all contigs
  --genomad-scores    geNomad aggregated classification of the hit contigs
  --genomad-plasmids  geNomad plasmid summary of the hit contigs (n_hallmarks, conjugation_genes)
  --mob-typer         mob_typer --multi on the hit contigs (optional)
  --paf               recruited reads -> their sample's contigs (minimap2 PAF), any number
  --taxa              hit read id -> Kraken2 taxid (HIT_READS; optional)
Writes <prefix>.tsv, <prefix>.fa and per_sample/<sample>.circles.fa (next to <prefix>; the Medaka input).
Hits of databases other than vfdb / plasmidfinder are listed as ARGs.
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

from extend_context import read_best_alignments

COLUMNS = ["contig", "sample", "length", "coverage", "reads", "plasmid_score", "chromosome_score", "hallmarks",
           "conjugation", "replicons", "mob_rep", "relaxase", "mpf", "mobility", "ARGs", "VFs", "read_taxa"]
ELEMENT = {"vfdb": "VFs", "plasmidfinder": "replicons"}


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def by_key(rows, key):
    return {r[key]: r for r in rows if r.get(key)}


def read_hits(path):
    """ABRicate tab output -> {contig: {"ARGs" | "VFs" | "replicons": set(gene)}}."""
    out = defaultdict(lambda: defaultdict(set))
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("#") or len(f) < 12:
                continue
            out[f[1]][ELEMENT.get(f[11].lower(), "ARGs")].add(f[5])
    return out


def read_fasta(path):
    seqs, name = {}, None
    if not path or not os.path.exists(path):
        return seqs
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def read_taxa(path):
    out = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1] not in ("", "0"):
                    out[f[0]] = f[1]
    return out


def clean(value):
    return "" if value in (None, "-") else value


def joined(values):
    return ";".join(sorted(v for v in values if v))


def candidates(info, hits, scores, plasmids, mob, best, taxa):
    reads, votes = Counter(), defaultdict(Counter)
    for read, aln in best.items():
        reads[aln[0]] += 1
        if read in taxa:
            votes[aln[0]][taxa[read]] += 1
    rows = []
    for r in info:
        c = r["contig"]
        if r["circular"] != "yes" or c not in hits:
            continue
        s, p, m, h = scores.get(c, {}), plasmids.get(c, {}), mob.get(c, {}), hits[c]
        conj = clean(p.get("conjugation_genes", ""))
        rows.append({"contig": c, "sample": r["sample"], "length": int(r["length"]), "coverage": r["coverage"],
                     "reads": reads[c], "plasmid_score": s.get("plasmid_score", ""),
                     "chromosome_score": s.get("chromosome_score", ""), "hallmarks": p.get("n_hallmarks", ""),
                     "conjugation": "" if conj == "NA" else conj,
                     "replicons": joined(h["replicons"]), "mob_rep": clean(m.get("rep_type(s)")),
                     "relaxase": clean(m.get("relaxase_type(s)")), "mpf": clean(m.get("mpf_type")),
                     "mobility": clean(m.get("predicted_mobility")), "ARGs": joined(h["ARGs"]), "VFs": joined(h["VFs"]),
                     "read_taxa": ",".join("%s:%d" % kv for kv in votes[c].most_common(3))})
    rows.sort(key=lambda x: (-x["length"], x["contig"]))
    return rows


def wrap(seq, width=80):
    return "".join(seq[i:i + width] + "\n" for i in range(0, len(seq), width))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--assembly-info", nargs="*", default=[])
    ap.add_argument("--hits", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--genomad-scores", required=True)
    ap.add_argument("--genomad-plasmids", required=True)
    ap.add_argument("--mob-typer", default="")
    ap.add_argument("--paf", nargs="*", default=[])
    ap.add_argument("--taxa", default="")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args(argv)
    info = [r for p in a.assembly_info for r in read_tsv(p)]
    rows = candidates(info, read_hits(a.hits), by_key(read_tsv(a.genomad_scores), "seq_name"),
                      by_key(read_tsv(a.genomad_plasmids), "seq_name"), by_key(read_tsv(a.mob_typer), "sample_id"),
                      read_best_alignments(a.paf), read_taxa(a.taxa))
    seqs = read_fasta(a.fasta)
    with open(a.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    per_sample = defaultdict(list)
    with open(a.out_prefix + ".fa", "w", newline="\n") as fo:
        for r in rows:
            seq = seqs.get(r["contig"], "")
            fo.write(">%s sample=%s length=%d replicons=%s ARGs=%s\n%s" % (
                r["contig"], r["sample"], len(seq), r["replicons"] or "-", r["ARGs"] or "-", wrap(seq)))
            per_sample[r["sample"]].append(">%s\n%s" % (r["contig"], wrap(seq)))
    out_dir = os.path.join(os.path.dirname(a.out_prefix), "per_sample")
    os.makedirs(out_dir, exist_ok=True)
    for sample, records in per_sample.items():
        with open(os.path.join(out_dir, sample + ".circles.fa"), "w", newline="\n") as fo:
            fo.write("".join(records))
    print("circular contigs: %d; with a hit: %d" % (sum(r["circular"] == "yes" for r in info), len(rows)),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
