#!/usr/bin/env python3
"""Per sample x panel species: Kraken2 clade reads, Bracken estimate and the next step.

  pathogen_candidates.py --panel panel_resolved.tsv --samples samples.tsv \
      --kreports <s>.kraken2.report.txt ... [--bracken <s>.tsv ...] --min-reads 10 \
      --out pathogen_candidates.tsv --samples-out pathogen_samples.tsv

decision (exactly one per sample x species):
  low_depth        sample flagged LOW_DEPTH by PREPROCESS: no reconstruction, no pooling
  detect           panel action = detect (reads + Bracken only)
  not_detected     no Kraken2 read in the species clade
  below_min_reads  1 .. min_reads-1 clade reads: reported LOW_COVERAGE, not extracted
  below_est_depth  >= min_reads clade reads, but the estimated depth (clade reads x the sample's mean
                   read length / reference genome size) is below --min-est-depth in the sample AND
                   summed over the patient's non-control samples: the coverage gate cannot pass, so
                   reported LOW_COVERAGE, not extracted
  extract          otherwise: extract (--include-children), gate, maybe assemble. Also when the depth
                   cannot be estimated (no seqkit stats / reference), as before the filter existed.
  control          a negative control that would be extracted: detection only (reported
                   DETECTED_ONLY) unless --gate-controls
est_depth = the estimate above (NA when unknown).
clade_reads = the clade count of the species line in the Kraken2 report (species + strains).
bracken_reads / bracken_fraction = Bracken new_est_reads / fraction_total_reads; 0 when the
species is absent from the sample's Bracken table, NA when the sample has no Bracken table.
genome_size = length of the panel reference genome (blank when not staged); the report's size ratio.
"""
import argparse
import csv
import gzip
import os
import sys

CAND_COLUMNS = ["sample", "patient", "is_control", "low_depth", "species", "slug", "tier", "action", "taxid",
                "clade_reads", "bracken_reads", "bracken_fraction", "est_depth", "decision", "genome_size"]
SAMPLE_COLUMNS = ["sample", "patient", "is_control", "low_depth", "n_reads", "classified_reads",
                  "unclassified_reads", "unclassified_pct"]


def sample_name(path, suffix):
    name = os.path.basename(str(path))
    return name[: -len(suffix)] if name.endswith(suffix) else os.path.splitext(name)[0]


def truthy(value):
    return str(value).strip().lower() in ("true", "1", "yes")


def read_kreport(path):
    """({taxid: clade reads}, unclassified, classified)"""
    clade, unclassified, classified = {}, 0, 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or not f[1].strip().isdigit():
                continue
            n, taxid = int(f[1]), f[4].strip()
            if f[3].strip() == "U":
                unclassified = n
                continue
            clade[taxid] = n
            if taxid == "1":
                classified = n
    return clade, unclassified, classified


def read_bracken(path):
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["taxonomy_id"]] = (row["new_est_reads"], row["fraction_total_reads"])
    return out


def decide(action, low_depth, clade, min_reads):
    if low_depth:
        return "low_depth"
    if action == "detect":
        return "detect"
    if clade == 0:
        return "not_detected"
    if clade < min_reads:
        return "below_min_reads"
    return "extract"


def genome_size(path):
    """Total sequence length of a (gzipped) FASTA."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return sum(len(line.strip()) for line in fh if not line.startswith(">"))


def read_mean_lengths(paths):
    """seqkit stats -T tables (<sample>.tsv): sample -> mean read length (sum_len / num_seqs)."""
    out = {}
    for path in paths:
        with open(path, newline="") as fh:
            row = next(csv.DictReader(fh, delimiter="	"), None)
        if row:
            n = int(str(row["num_seqs"]).replace(",", ""))
            bases = int(str(row["sum_len"]).replace(",", ""))
            if n > 0:
                out[sample_name(path, ".tsv")] = bases / n
    return out


def run(panel_path, samples_path, kreports, brackens, min_reads, out_path, samples_out,
        seqkit_stats=(), refdir=None, min_est_depth=0.0, gate_controls=False):
    with open(panel_path, newline="") as fh:
        panel = list(csv.DictReader(fh, delimiter="\t"))
    with open(samples_path, newline="") as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))
    reports = {sample_name(p, ".kraken2.report.txt"): p for p in kreports}
    bracken = {sample_name(p, ".tsv"): read_bracken(p) for p in brackens}
    mean_len = read_mean_lengths(seqkit_stats)
    sizes = {}
    for p in panel:
        ref = os.path.join(refdir, "%s.fna.gz" % p["reference_accession"]) if refdir and p["reference_accession"] else None
        if ref and os.path.isfile(ref):
            sizes[p["species"]] = genome_size(ref)
    cand_rows, sample_rows = [], []
    for s in samples:
        sid = s["sample"]
        if sid not in reports:
            raise ValueError("no Kraken2 report for sample %s" % sid)
        clade, unclassified, classified = read_kreport(reports[sid])
        total = unclassified + classified
        low = truthy(s["low_depth"])
        sample_rows.append({
            "sample": sid, "patient": s["patient"], "is_control": s["is_control"], "low_depth": s["low_depth"],
            "n_reads": s["n_reads"], "classified_reads": classified, "unclassified_reads": unclassified,
            "unclassified_pct": "%.2f" % (100.0 * unclassified / total if total else 0.0),
        })
        for p in panel:
            n = clade.get(p["taxid"], 0)
            if sid in bracken:
                b_reads, b_frac = bracken[sid].get(p["taxid"], ("0", "0"))
            else:
                b_reads, b_frac = "NA", "NA"
            g, m = sizes.get(p["species"]), mean_len.get(sid)
            est = n * m / g if g and m else None
            decision = decide(p["action"], low, n, min_reads)
            if decision == "extract" and truthy(s["is_control"]) and not gate_controls:
                decision = "control"   # negative controls: detection only unless --gate-controls
            cand_rows.append({
                "sample": sid, "patient": s["patient"], "is_control": s["is_control"], "low_depth": s["low_depth"],
                "species": p["species"], "slug": p["slug"], "tier": p["tier"], "action": p["action"],
                "taxid": p["taxid"], "clade_reads": n, "bracken_reads": b_reads, "bracken_fraction": b_frac,
                "est_depth": "NA" if est is None else "%.2f" % est,
                "decision": decision, "genome_size": g or "",
            })
    apply_depth_filter(cand_rows, min_est_depth)
    for path, columns, rows in ((out_path, CAND_COLUMNS, cand_rows), (samples_out, SAMPLE_COLUMNS, sample_rows)):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    return cand_rows


def apply_depth_filter(rows, min_est_depth):
    """extract -> below_est_depth when neither the sample nor its patient's pooled samples can reach
    min_est_depth. Pooled = the patient's non-control samples that are extracted for that species."""
    if not min_est_depth or min_est_depth <= 0:
        return
    pooled = {}
    for r in rows:
        if r["decision"] == "extract" and r["est_depth"] != "NA" and r["patient"] and not truthy(r["is_control"]):
            key = (r["patient"], r["species"])
            pooled[key] = pooled.get(key, 0.0) + float(r["est_depth"])
    for r in rows:
        if r["decision"] != "extract" or r["est_depth"] == "NA":
            continue
        if float(r["est_depth"]) >= min_est_depth:
            continue
        if r["patient"] and not truthy(r["is_control"]) and pooled.get((r["patient"], r["species"]), 0.0) >= min_est_depth:
            continue
        r["decision"] = "below_est_depth"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--kreports", nargs="+", required=True)
    ap.add_argument("--bracken", nargs="*", default=[])
    ap.add_argument("--min-reads", type=int, default=10)
    ap.add_argument("--seqkit-stats", nargs="*", default=[], help="seqkit stats -T tables (<sample>.tsv)")
    ap.add_argument("--refdir", default=None, help="panel reference genomes (<accession>.fna.gz)")
    ap.add_argument("--min-est-depth", type=float, default=0.0,
                    help="skip extraction below this estimated depth (sample and patient-pooled); 0 = off")
    ap.add_argument("--gate-controls", action="store_true",
                    help="extract and gate negative controls like samples (default: detection only)")
    ap.add_argument("--out", default="pathogen_candidates.tsv")
    ap.add_argument("--samples-out", default="pathogen_samples.tsv")
    args = ap.parse_args(argv)
    try:
        rows = run(args.panel, args.samples, args.kreports, args.bracken, args.min_reads, args.out, args.samples_out,
                   args.seqkit_stats, args.refdir, args.min_est_depth, args.gate_controls)
    except (ValueError, KeyError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    n_extract = sum(r["decision"] == "extract" for r in rows)
    n_skip = sum(r["decision"] == "below_est_depth" for r in rows)
    print("%d sample x species rows, %d to extract, %d below --min-est-depth" % (len(rows), n_extract, n_skip),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
