#!/usr/bin/env python3
"""Species confirmation of one assembly from skani dist output.

  species_confirm.py --unit <id> --species "<panel species>" --ani-min 95 --skani <id>.skani.tsv \
      --references panel_references.tsv [--confusable-delta 0.5] --out <id>.species.tsv

skani compared the assembly with the panel reference and its confusables (Ref_file =
<accession>.fna.gz). The best hit (highest ANI, ties -> higher query aligned fraction) must be the
panel species with ANI >= ani-min, otherwise SPECIES_NOT_CONFIRMED and no typing is reported
under the panel species name.

reason: confirmed | ani_below_min | best_hit_other_species | unresolved_confusable | no_skani_hit.
unresolved_confusable: the best hit is a confusable relative, but the panel species' own reference
is within --confusable-delta ANI points of it and itself >= ani-min (e.g. E. coli vs Shigella):
ANI cannot separate them, marker genes / MLST are needed. panel_ani = the panel reference's ANI.
"""
import argparse
import csv
import os
import sys

OUT_COLUMNS = ["unit_id", "panel_species", "best_species", "best_accession", "ani", "af_query", "confirmed", "reason",
               "panel_ani"]
SUFFIXES = (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta")


def accession_of(ref_file):
    name = os.path.basename(ref_file)
    for suffix in SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def confirm(unit, species, ani_min, skani, references, delta=0.0):
    with open(references, newline="") as fh:
        by_acc = {r["accession"].strip(): r["species"].strip() for r in csv.DictReader(fh, delimiter="\t")}
    with open(skani, newline="") as fh:
        hits = list(csv.DictReader(fh, delimiter="\t"))
    row = {"unit_id": unit, "panel_species": species, "best_species": "", "best_accession": "",
           "ani": "", "af_query": "", "confirmed": "false", "reason": "no_skani_hit", "panel_ani": ""}
    if not hits:
        return row
    best = max(hits, key=lambda h: (float(h["ANI"]), float(h["Align_fraction_query"])))
    acc = accession_of(best["Ref_file"])
    row.update(best_species=by_acc.get(acc, acc), best_accession=acc, ani=best["ANI"], af_query=best["Align_fraction_query"])
    own = [h for h in hits if by_acc.get(accession_of(h["Ref_file"])) == species]
    own_ani = max((float(h["ANI"]) for h in own), default=None)
    if own_ani is not None:
        row["panel_ani"] = "%.2f" % own_ani
    if row["best_species"] != species:
        close = own_ani is not None and own_ani >= ani_min and float(best["ANI"]) - own_ani <= delta
        row["reason"] = "unresolved_confusable" if close else "best_hit_other_species"
    elif float(best["ANI"]) < ani_min:
        row["reason"] = "ani_below_min"
    else:
        row.update(confirmed="true", reason="confirmed")
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--unit", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--ani-min", type=float, required=True)
    ap.add_argument("--skani", required=True)
    ap.add_argument("--references", required=True)
    ap.add_argument("--confusable-delta", type=float, default=0.0,
                    help="ANI points within which a confusable best hit leaves the species unresolved")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    row = confirm(args.unit, args.species, args.ani_min, args.skani, args.references, args.confusable_delta)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
