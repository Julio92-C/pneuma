#!/usr/bin/env python3
"""Genomic context of ARG hits (or VF / MGE hits, --abricate-db): plasmid / chromosomal_mge / chromosome / unresolved.

Inputs are pooled over samples: every sequence id is "<sample>__<id>" (an ONT read for --route reads, an
assembled contig for --route contigs). One output row per ABRicate CARD hit, classified from geNomad
(aggregated scores, proviruses, gene annotations) and MOB-suite. MOB-suite only adds plasmid evidence
(mob_recon molecule_type, or a replicon found by mob_typer on the same read): its "chromosome" label
means "no plasmid match", not a chromosomal call.
"""
import argparse
import csv
import os
import re
import sys

CATEGORIES = ["plasmid", "chromosomal_mge", "chromosome", "unresolved"]
LOCUS_COLUMNS = ["sample", "route", "seq", "seq_len", "gene", "drug_class", "start", "end", "category",
                 "reason", "plasmid_score", "chromosome_score", "virus_score", "in_provirus", "mge_genes",
                 "mob_molecule", "mob_detail"]
MGE_RE = re.compile(r"integrase|transposase|recombinase|resolvase|invertase|relaxase|mobili[sz]ation|"
                    r"conjugal|conjugative|excisionase|insertion sequence", re.I)
GENE_ID_RE = re.compile(r"^(.+)_(\d+)$")
PROVIRUS_RE = re.compile(r"^(.+)\|provirus_(\d+)_(\d+)$")
WHOLE = (1, 10 ** 12)   # a sequence that is a virus as a whole


def exists(path):
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def read_tsv(path):
    if not exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sample_of(seq):
    if "__" not in seq:
        raise ValueError("sequence id without a '<sample>__' prefix: %s" % seq)
    return seq.split("__", 1)[0]


def read_abricate(path, db="card"):
    """ABRicate tab or --csv output (headers may repeat) -> hits of one database."""
    if not exists(path):
        return []
    hits, header = [], None
    with open(path, newline="") as fh:
        for line in fh:
            if not line.strip():
                continue
            sep = "\t" if "\t" in line else ","
            row = next(csv.reader([line.rstrip("\r\n")], delimiter=sep))
            if row[0].startswith("#FILE"):
                header = row
                continue
            if header is None:
                raise ValueError("ABRicate file without a #FILE header: %s" % path)
            r = dict(zip(header, row))
            if r.get("DATABASE", "").lower() != db:
                continue
            a, b = int(r["START"]), int(r["END"])
            hits.append({"seq": r["SEQUENCE"], "start": min(a, b), "end": max(a, b), "gene": r["GENE"],
                         "drug_class": r.get("RESISTANCE", "")})
    return hits


def read_lengths(path):
    """seqkit fx2tab -n -i -l output: id<TAB>length[<TAB>]."""
    out = {}
    if exists(path):
        for line in open(path):
            vals = [v for v in line.rstrip("\r\n").split("\t") if v.strip()]
            if len(vals) >= 2:
                out[vals[0].strip()] = int(vals[-1])
    return out


def read_genomad_scores(path):
    """aggregated_classification.tsv -> {seq: (plasmid, chromosome, virus)}."""
    return {r["seq_name"]: (float(r["plasmid_score"]), float(r["chromosome_score"]), float(r["virus_score"]))
            for r in read_tsv(path)}


def read_proviruses(path):
    """virus_summary.tsv -> {seq: [(start, end), ...]}."""
    out = {}
    for r in read_tsv(path):
        m = PROVIRUS_RE.match(r["seq_name"])
        seq, span = (m.group(1), (int(m.group(2)), int(m.group(3)))) if m else (r["seq_name"], WHOLE)
        out.setdefault(seq, []).append(span)
    return out


def read_mge_genes(path):
    """genes.tsv -> {seq: [(start, end, label)]} for genes with a CONJscan hit or a mobile-element description."""
    out = {}
    for r in read_tsv(path):
        conj = r.get("annotation_conjscan", "") or ""
        conj = "" if conj == "NA" else conj
        desc = r.get("annotation_description", "") or ""
        if not conj and not MGE_RE.search(desc):
            continue
        m = GENE_ID_RE.match(r["gene"])
        if m:
            out.setdefault(m.group(1), []).append((int(r["start"]), int(r["end"]), conj or desc))
    return out


def read_mob_contigs(path):
    """mob_recon contig_report.txt (concatenated over samples) -> {contig: (molecule_type, primary_cluster_id)}."""
    out = {}
    for r in read_tsv(path):
        if r.get("contig_id") in (None, "contig_id"):
            continue
        cluster = r.get("primary_cluster_id", "")
        out[r["contig_id"].split()[0]] = (r.get("molecule_type", ""), "" if cluster == "-" else cluster)
    return out


def read_mob_typer(path):
    """mob_typer --multi results -> {seq: ("plasmid" if a replicon is on the sequence else "", replicons)}."""
    out = {}
    for r in read_tsv(path):
        if r.get("sample_id") in (None, "sample_id"):
            continue
        rep = r.get("rep_type(s)", "-") or "-"
        out[r["sample_id"].split()[0]] = ("plasmid" if rep != "-" else "", "" if rep == "-" else rep)
    return out


def read_plasmid_evidence(path):
    """plasmid_summary.tsv -> sequences with a geNomad plasmid hallmark or conjugation gene."""
    return {r["seq_name"] for r in read_tsv(path)
            if (r.get("n_hallmarks") or "0") not in ("0", "NA", "") or (r.get("conjugation_genes") or "NA") != "NA"}


def classify(scores, mob_molecule="", in_provirus=False, mge_near=False, flank_ok=True,
             plasmid_min=0.7, chrom_min=0.7, plasmid_supported=True):
    """scores = (plasmid, chromosome, virus) from geNomad, or None when geNomad has no row.

    plasmid_supported=False: the sequence has no plasmid hallmark / conjugation gene, so a geNomad plasmid
    score alone does not make a plasmid call (a replicon from MOB-suite / PlasmidFinder still does)."""
    if not flank_ok:
        return "unresolved", "short_flank"
    pl, ch, _vi = scores if scores else (0.0, 0.0, 0.0)
    g_pl, g_ch = pl >= plasmid_min, ch >= chrom_min
    mob_pl = mob_molecule == "plasmid"
    score_only = g_pl and not plasmid_supported and not mob_pl
    g_pl = g_pl and not score_only
    if mob_pl and g_ch:
        return "unresolved", "conflict_mob_plasmid_genomad_chromosome"
    if g_pl or mob_pl:
        return "plasmid", "+".join(n for n, v in (("genomad", g_pl), ("mobsuite", mob_pl)) if v)
    if in_provirus:
        return "chromosomal_mge", "provirus"
    if g_ch and mge_near:
        return "chromosomal_mge", "mge_genes_nearby"
    if g_ch:
        return "chromosome", "genomad"
    if score_only:
        return "unresolved", "plasmid_score_only"
    return "unresolved", "low_scores" if scores else "no_scores"


def overlaps(a, b, pad=0):
    return a[0] - pad <= b[1] and b[0] <= a[1] + pad


def loci(route, hits, lengths, scores, proviruses, mge, mob, min_flank=1000, mge_window=10000,
         plasmid_min=0.7, chrom_min=0.7, plasmid_evidence=None):
    """plasmid_evidence: None = geNomad plasmid scores are used as they are; a set = only those sequences may get
    a plasmid call from the geNomad score (see classify)."""
    rows = []
    for h in hits:
        seq, span = h["seq"], (h["start"], h["end"])
        n = lengths.get(seq)
        if n is None:
            raise ValueError("no length for %s" % seq)
        flank_ok = span[0] - 1 >= min_flank and n - span[1] >= min_flank
        s = scores.get(seq)
        in_pro = any(overlaps(span, p) for p in proviruses.get(seq, []))
        near = sorted({lab for a, b, lab in mge.get(seq, [])
                       if overlaps(span, (a, b), mge_window) and not overlaps(span, (a, b))})
        molecule, detail = mob.get(seq, ("", ""))
        supported = plasmid_evidence is None or seq in plasmid_evidence
        cat, reason = classify(s, molecule, in_pro, bool(near), flank_ok, plasmid_min, chrom_min, supported)
        rows.append({"sample": sample_of(seq), "route": route, "seq": seq, "seq_len": n, "gene": h["gene"],
                     "drug_class": h["drug_class"], "start": span[0], "end": span[1], "category": cat,
                     "reason": reason,
                     "plasmid_score": "%.3f" % s[0] if s else "",
                     "chromosome_score": "%.3f" % s[1] if s else "",
                     "virus_score": "%.3f" % s[2] if s else "",
                     "in_provirus": "yes" if in_pro else "no", "mge_genes": ";".join(near),
                     "mob_molecule": molecule, "mob_detail": detail})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--route", required=True, choices=["reads", "contigs"])
    ap.add_argument("--abricate", required=True, help="ABRicate output over all samples")
    ap.add_argument("--abricate-db", default="card", help="ABRicate database whose hits are classified (card, vfdb, plasmidfinder...)")
    ap.add_argument("--lengths", required=True, help="seqkit fx2tab -n -i -l of the classified sequences")
    ap.add_argument("--genomad-scores", default="", help="<prefix>_aggregated_classification.tsv")
    ap.add_argument("--genomad-virus", default="", help="<prefix>_virus_summary.tsv")
    ap.add_argument("--genomad-genes", default="", help="<prefix>_genes.tsv (annotate step)")
    ap.add_argument("--mob-contigs", default="", help="mob_recon contig_report.txt (route contigs)")
    ap.add_argument("--mob-typer", default="", help="mob_typer --multi results (route reads)")
    ap.add_argument("--genomad-plasmids", default="",
                    help="<prefix>_plasmid_summary.tsv: when given, a geNomad plasmid call needs a plasmid hallmark or "
                         "conjugation gene on the sequence (or a replicon from MOB-suite / PlasmidFinder)")
    ap.add_argument("--min-flank", type=int, default=1000, help="bp required on both sides of the gene")
    ap.add_argument("--mge-window", type=int, default=10000, help="bp around the gene searched for MGE genes")
    ap.add_argument("--plasmid-min", type=float, default=0.7)
    ap.add_argument("--chrom-min", type=float, default=0.7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    mob = read_mob_contigs(a.mob_contigs) if a.mob_contigs else read_mob_typer(a.mob_typer)
    rows = loci(a.route, read_abricate(a.abricate, a.abricate_db.lower()), read_lengths(a.lengths), read_genomad_scores(a.genomad_scores),
                read_proviruses(a.genomad_virus), read_mge_genes(a.genomad_genes), mob,
                a.min_flank, a.mge_window, a.plasmid_min, a.chrom_min,
                read_plasmid_evidence(a.genomad_plasmids) if a.genomad_plasmids else None)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=LOCUS_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print("%d ARG loci (%s): %s" % (len(rows), a.route, ", ".join(
        "%s=%d" % (c, sum(r["category"] == c for r in rows)) for c in CATEGORIES)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
