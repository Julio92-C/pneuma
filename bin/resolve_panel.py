#!/usr/bin/env python3
"""Resolve the PNEUMA pathogen panel (species names) to NCBI taxids.

  resolve_panel.py --panel pathogen_panel.csv --references panel_references.tsv \
                   --taxdump <dir with nodes.dmp, names.dmp> --out panel_resolved.tsv

Panel columns: species, tier, who_category, action, reference_accession, confusable_with, ani_min,
amrfinder_organism (optional). `confusable_with` is a ';'-separated list of species names, each
needing an accession in panel_references.tsv (species, accession, role).

Names match NCBI scientific names first (case-insensitive), then synonyms / equivalent names, so a
renamed species keeps its taxid ('Haemophilus parainfluenzae' -> 729 'Parahaemophilus
parainfluenzae', with a warning). Every problem is reported at once; unresolved names get the
closest species names as suggestions; exit 1. `taxids` = the species taxid plus all descendants
(what extract_kraken_reads.py --include-children collects).
"""
import argparse
import csv
import difflib
import sys
from collections import defaultdict

TIERS = ("1", "2", "3", "D", "R")
ACTIONS = ("reconstruct", "detect")
NAME_CLASSES = ("scientific name", "synonym", "equivalent name", "genbank synonym", "basionym")
PANEL_COLUMNS = ["species", "tier", "who_category", "action", "reference_accession", "confusable_with", "ani_min"]
OUT_COLUMNS = ["species", "slug", "tier", "who_category", "action", "taxid", "scientific_name",
               "matched_name_class", "taxids", "reference_accession", "confusable_with",
               "confusable_accessions", "ani_min", "amrfinder_organism"]
DEFAULT_ANI_MIN = 95.0


def norm(name):
    return " ".join(name.replace("_", " ").split()).lower()


def slug(name):
    return "_".join("".join(c if c.isalnum() else " " for c in name).split())


def dmp(line):
    return [f.strip() for f in line.rstrip("\n").split("|")]


def read_nodes(path):
    parent, rank = {}, {}
    with open(path) as fh:
        for line in fh:
            f = dmp(line)
            if len(f) >= 3 and f[0]:
                parent[f[0]] = f[1]
                rank[f[0]] = f[2]
    return parent, rank


def read_names(path, wanted):
    """{normalised wanted name: [(taxid, class)]} and {taxid: scientific name} (all taxa)."""
    hits, sci = defaultdict(list), {}
    with open(path) as fh:
        for line in fh:
            f = dmp(line)
            if len(f) < 4:
                continue
            if f[3] == "scientific name":
                sci[f[0]] = f[1]
            key = norm(f[1])
            if key in wanted and f[3] in NAME_CLASSES:
                hits[key].append((f[0], f[3]))
    return hits, sci


def build_children(parent):
    children = defaultdict(list)
    for child, par in parent.items():
        if child != par:
            children[par].append(child)
    return children


def descendants(taxid, children):
    out, stack = [], [taxid]
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(children.get(node, []))
    return sorted(out, key=int)


def resolve(name, hits, rank):
    """(taxid, name class, None) or (None, None, problem)."""
    cands = hits.get(norm(name), [])
    sci = sorted({t for t, c in cands if c == "scientific name"}, key=int)
    pool = sci or sorted({t for t, _c in cands}, key=int)
    if not pool:
        return None, None, "not found in the taxonomy"
    if len(pool) > 1:
        return None, None, "ambiguous, matches taxids " + ", ".join(pool)
    taxid = pool[0]
    if rank.get(taxid) != "species":
        return None, None, "taxid %s has rank '%s', expected species" % (taxid, rank.get(taxid))
    cls = "scientific name" if sci else sorted(c for t, c in cands if t == taxid)[0]
    return taxid, cls, None


def suggestions(name, sci, rank):
    head = norm(name)[:2]
    pool = sorted({n for t, n in sci.items() if rank.get(t) == "species" and n[:2].lower() == head})
    return difflib.get_close_matches(name, pool, n=3, cutoff=0.6)


def read_references(path):
    with open(path, newline="") as fh:
        return {r["species"].strip(): r["accession"].strip() for r in csv.DictReader(fh, delimiter="\t")}


def run(panel_path, references_path, taxdump, out_path, log=sys.stderr):
    with open(panel_path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in PANEL_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError("pathogen panel %s lacks column(s): %s" % (panel_path, ", ".join(missing)))
        rows = [{k: (v or "").strip() for k, v in r.items() if k} for r in reader]
    refs = read_references(references_path)
    problems, seen = [], set()
    for i, r in enumerate(rows, 2):
        where = "%s line %d (%s)" % (panel_path, i, r["species"] or "<empty species>")
        if not r["species"]:
            problems.append("%s: empty species" % where)
        if r["species"] in seen:
            problems.append("%s: duplicated species" % where)
        seen.add(r["species"])
        if r["tier"] not in TIERS:
            problems.append("%s: tier '%s' not in %s" % (where, r["tier"], "/".join(TIERS)))
        if r["action"] not in ACTIONS:
            problems.append("%s: action '%s' not in %s" % (where, r["action"], "/".join(ACTIONS)))
        if r["action"] == "reconstruct":
            if not r["reference_accession"]:
                problems.append("%s: reconstruct needs reference_accession" % where)
            elif refs.get(r["species"]) != r["reference_accession"]:
                problems.append("%s: reference_accession %s does not match %s (%s)" % (
                    where, r["reference_accession"], references_path, refs.get(r["species"], "species absent")))
            try:
                r["ani_min"] = str(float(r["ani_min"])) if r["ani_min"] else str(DEFAULT_ANI_MIN)
            except ValueError:
                problems.append("%s: ani_min '%s' is not a number" % (where, r["ani_min"]))
        conf = [c.strip() for c in r.get("confusable_with", "").split(";") if c.strip()]
        absent = [c for c in conf if c not in refs]
        if absent:
            problems.append("%s: confusable species without an accession in %s: %s"
                            % (where, references_path, ", ".join(absent)))
        r["_confusable_accessions"] = [refs[c] for c in conf if c in refs]
    parent, rank = read_nodes("%s/nodes.dmp" % taxdump)
    hits, sci = read_names("%s/names.dmp" % taxdump, {norm(r["species"]) for r in rows})
    children = build_children(parent)
    out = []
    for r in rows:
        taxid, cls, problem = resolve(r["species"], hits, rank)
        if problem:
            sug = suggestions(r["species"], sci, rank)
            problems.append("species '%s': %s%s" % (r["species"], problem,
                                                    ("; did you mean: " + ", ".join(sug)) if sug else ""))
            continue
        if cls != "scientific name":
            log.write("WARNING: '%s' matched the %s of taxid %s; the current NCBI name is '%s'\n"
                      % (r["species"], cls, taxid, sci.get(taxid, "?")))
        out.append({
            "species": r["species"], "slug": slug(r["species"]), "tier": r["tier"],
            "who_category": r["who_category"], "action": r["action"], "taxid": taxid,
            "scientific_name": sci.get(taxid, ""), "matched_name_class": cls,
            "taxids": ";".join(descendants(taxid, children)),
            "reference_accession": r["reference_accession"], "confusable_with": r.get("confusable_with", ""),
            "confusable_accessions": ";".join(r["_confusable_accessions"]),
            "ani_min": r["ani_min"] if r["action"] == "reconstruct" else "",
            "amrfinder_organism": r.get("amrfinder_organism", ""),
        })
    if problems:
        raise ValueError("Cannot resolve the pathogen panel:\n  - " + "\n  - ".join(problems))
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--references", required=True)
    ap.add_argument("--taxdump", required=True)
    ap.add_argument("--out", default="panel_resolved.tsv")
    args = ap.parse_args(argv)
    try:
        rows = run(args.panel, args.references, args.taxdump, args.out)
    except ValueError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    print("resolved %d panel species" % len(rows), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
