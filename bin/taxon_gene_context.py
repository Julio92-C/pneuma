#!/usr/bin/env python3
"""Every read hit of every context database with its genomic context and its read's Kraken2 taxon (MERIDIAN hand-off).

Inputs
  --hit-context  hit_context_<db>.tsv (extend_context.py, local assembly), any number
  --read-loci    loci_reads_<db>.tsv (arg_context.py, evidence rule): the read calls of a database without a
                 hit_context table (--context_local_assembly false); their contig columns are empty (no_contig_match)
  --taxa         hit read id -> Kraken2 taxid (HIT_READS, pooled; optional)
  --taxdump      directory with nodes.dmp and names.dmp
Taxon fields of the hit read: kraken_taxid, kraken_rank, kraken_name, species, genus, family, phylum and taxon_level
(the most specific of species / genus / family / above_family; unclassified for taxid 0 or none, and for a taxid
that is not in the taxdump, e.g. a Kraken2 database and a taxdump from different releases).
Host vote of a contig: each hit read placed on the contig votes once, for its species (its genus when it has no
species); contig_taxon is the most-voted name, contig_vote its share of the contig's hit reads and contig_vote_reads
the reads behind it. weak_vote = yes when that share is below --min-vote (0.5), fewer than --min-vote-reads (2)
reads voted for it, or another taxon has as many votes (a tie; contig_taxon is then the first name alphabetically); NA for a hit that is not on a contig (the read's own taxon is the only evidence).
Writes --out (hit_context.tsv) and --resolution (taxon_resolution.tsv: per element and final context, the share of
hits whose read is named to species, genus, family, only above family, or unclassified).
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bin_recovery import ancestor_at, read_taxonomy  # noqa: E402
from extend_context import COLUMNS as CONTEXT_COLUMNS, extend, read_tsv  # noqa: E402

ELEMENT = {"card": "AMR", "vfdb": "VFs", "plasmidfinder": "MGEs"}
TAXON_COLUMNS = ["kraken_taxid", "kraken_rank", "kraken_name", "species", "genus", "family", "phylum", "taxon_level",
                 "contig_taxon", "contig_vote", "contig_vote_reads", "weak_vote"]
COLUMNS = CONTEXT_COLUMNS[:3] + ["element"] + CONTEXT_COLUMNS[3:] + TAXON_COLUMNS
LEVELS = ["species", "genus", "family", "above_family", "unclassified"]
RESOLUTION_COLUMNS = ["element", "context", "hits"] + ["%s_pct" % lv for lv in LEVELS]
HIT_CONTEXT = re.compile(r"hit_context_(.+)\.tsv$")
READ_LOCI = re.compile(r"loci_reads_(.+)\.tsv$")


def db_of(path, pattern):
    name = os.path.basename(path)
    if name.endswith("_score_only.tsv"):
        raise ValueError("%s: a score-only table; pass the evidence-rule loci_reads_<db>.tsv" % path)
    m = pattern.search(name)
    if not m:
        raise ValueError("%s: expected a file name matching %s" % (path, pattern.pattern))
    return m.group(1).lower()


def load(hit_context_paths, read_loci_paths):
    """Rows of the hit_context tables, then read-only rows for each database that has no hit_context table."""
    rows, done = [], set()
    for p in hit_context_paths:
        done.add(db_of(p, HIT_CONTEXT))
        rows.extend(read_tsv(p))
    for p in read_loci_paths:
        db = db_of(p, READ_LOCI)
        if db not in done:
            rows.extend(extend(read_tsv(p), [], {}, db))
    return rows


def read_taxa(path):
    out = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    out[f[0]] = f[1]
    return out


def lineage(taxid, taxonomy):
    """Kraken2 taxid -> taxon fields (names empty when the taxdump does not know the taxid)."""
    parent, rank, name = taxonomy
    fields = {"kraken_taxid": taxid, "kraken_rank": rank.get(taxid, ""), "kraken_name": name.get(taxid, "")}
    for level in ("species", "genus", "family", "phylum"):
        anc = ancestor_at(taxid, level, parent, rank) if taxid else None
        fields[level] = name.get(anc, "") if anc else ""
    if taxid in ("", "0") or taxid not in parent:
        fields["taxon_level"] = "unclassified"
    else:
        fields["taxon_level"] = next((lv for lv in ("species", "genus", "family") if fields[lv]), "above_family")
    return fields


def host_votes(rows, fields):
    """contig -> (taxon, share, reads) from the distinct hit reads placed on it."""
    voters = defaultdict(dict)
    for r in rows:
        if r.get("contig"):
            t = fields[r["seq"]]
            voters[r["contig"]][r["seq"]] = t["species"] or t["genus"]
    out = {}
    for contig, votes in voters.items():
        counts = Counter(v for v in votes.values() if v)
        if not counts:
            out[contig] = ("", 0.0, 0, False)
            continue
        ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
        taxon, n = ranked[0]
        out[contig] = (taxon, n / len(votes), n, len(ranked) > 1 and ranked[1][1] == n)
    return out


def annotate(rows, taxa, taxonomy, min_vote=0.5, min_vote_reads=2):
    fields = {}
    for r in rows:
        if r["seq"] not in fields:
            fields[r["seq"]] = lineage(taxa.get(r["seq"], ""), taxonomy)
    votes = host_votes(rows, fields)
    out = []
    for r in rows:
        row = dict(r)
        row["element"] = ELEMENT.get(r["database"], r["database"])
        row.update(fields[r["seq"]])
        if r.get("contig"):
            taxon, share, n, tied = votes[r["contig"]]
            row.update(contig_taxon=taxon, contig_vote="%.3f" % share, contig_vote_reads=n,
                       weak_vote="yes" if tied or share < min_vote or n < min_vote_reads else "no")
        else:
            row.update(contig_taxon="", contig_vote="", contig_vote_reads="", weak_vote="NA")
        out.append(row)
    return out


def resolution(rows):
    acc = defaultdict(Counter)
    for r in rows:
        for ctx in (r["final_context"], "all"):
            acc[(r["element"], ctx)][r["taxon_level"]] += 1
    out = []
    for element, ctx in sorted(acc):
        c = acc[(element, ctx)]
        n = sum(c.values())
        row = {"element": element, "context": ctx, "hits": n}
        row.update(("%s_pct" % lv, "%.1f" % (100.0 * c[lv] / n)) for lv in LEVELS)
        out.append(row)
    return out


def write(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore",
                           restval="")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hit-context", nargs="*", default=[], help="hit_context_<db>.tsv tables")
    ap.add_argument("--read-loci", nargs="*", default=[], help="loci_reads_<db>.tsv tables (layer 1 only)")
    ap.add_argument("--taxa", default="", help="hit read id <TAB> Kraken2 taxid")
    ap.add_argument("--taxdump", required=True, help="directory with nodes.dmp and names.dmp")
    ap.add_argument("--min-vote", type=float, default=0.5, help="host vote share below which weak_vote = yes")
    ap.add_argument("--min-vote-reads", type=int, default=2, help="voting reads below which weak_vote = yes")
    ap.add_argument("--out", required=True)
    ap.add_argument("--resolution", required=True)
    a = ap.parse_args(argv)
    rows = annotate(load(a.hit_context, a.read_loci), read_taxa(a.taxa), read_taxonomy(a.taxdump), a.min_vote,
                    a.min_vote_reads)
    write(a.out, COLUMNS, rows)
    write(a.resolution, RESOLUTION_COLUMNS, resolution(rows))
    unknown = set(r["kraken_taxid"] for r in rows if r["kraken_taxid"] not in ("", "0") and not r["kraken_name"]
                  and r["taxon_level"] == "unclassified")
    if unknown:
        print("WARNING: %d read taxids are not in the taxdump; counted as unclassified" % len(unknown),
              file=sys.stderr)
    on = [r for r in rows if r["weak_vote"] != "NA"]
    print("%d hits; %d on a contig, %d of them with a weak host vote" % (
        len(rows), len(on), sum(r["weak_vote"] == "yes" for r in on)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
