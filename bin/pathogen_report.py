#!/usr/bin/env python3
"""PNEUMA pathogen report: one primary status per sample x panel species (spec 5.6), plus patient-level
rows for pooled reconstructions, a cohort TSV and one HTML page per sample.

  pathogen_report.py --panel panel_resolved.tsv --candidates pathogen_candidates.tsv \
      --samples pathogen_samples.tsv --units units.tsv --results <dir> --panel-file pathogen_panel.csv \
      [--manifest database_manifest.tsv] --pipeline-version 0.3.0 \
      --gate-min-depth 10 --gate-min-breadth 0.5 --outdir .

<dir>/ holds the per-unit results (a missing sub-directory = nothing of that kind was produced):
  gate/<unit>.gate.tsv  species/<unit>.species.tsv  assemblies/<unit>.fasta  mlst/<unit>.tsv
  amrfinderplus/<unit>.tsv  abricate/<unit>.{vfdb,plasmidfinder}.txt  mobsuite/<unit>.mobtyper_results.txt
  mixed/<unit>.mixed.tsv  checkm2/*_checkm2_report.tsv  quast/*_quast.tsv  functional/<sample>_abricate.csv
A missing gate / species / mixed / mlst result for a unit that needs one is a wiring error (exit 1).
A pooled (patient-level) result is shown under each contributing sample as
"<STATUS> (pooled, patient <id>)" and in a patient section; it never replaces the sample's status.
"""
import argparse
import csv
import hashlib
import html
import os
import re
import sys
from collections import defaultdict

STATUSES = ["NOT_DETECTED", "DETECTED_ONLY", "LOW_DEPTH", "LOW_COVERAGE", "ASSEMBLY_FAILED",
            "SPECIES_NOT_CONFIRMED", "MIXED_STRAINS", "TYPED", "PARTIAL", "NO_SCHEME"]
COLUMNS = ["level", "sample", "patient", "contributing_samples", "species", "tier", "who_category", "action", "taxid",
           "status", "pooled_result", "reads", "bracken_reads", "bracken_fraction", "depth", "breadth", "contigs",
           "total_length", "n50", "completeness", "contamination", "ani", "ani_reference", "scheme", "st", "alleles",
           "amr_genes", "vf_genes", "mge_genes", "gene_source", "mixed_het_per_kb", "notes"]
SECTIONS = [("Tier 1-3 priority pathogens", ("1", "2", "3")), ("Detect-only (Tier D)", ("D",)),
            ("Tier R commensal AMR/VF reservoirs", ("R",))]
STATUS_HELP = [
    ("TYPED", "assembled, species confirmed (skani), complete MLST profile"),
    ("PARTIAL", "assembled and confirmed; MLST loci missing, novel or ambiguous"),
    ("NO_SCHEME", "assembled and confirmed; no PubMLST scheme for this species"),
    ("MIXED_STRAINS", "pooled assembly with minor-allele sites above the threshold: more than one strain, no ST"),
    ("SPECIES_NOT_CONFIRMED", "best ANI hit is another species, or ANI below ani_min: no typing under this name"),
    ("ASSEMBLY_FAILED", "passed the coverage gate but Flye/Medaka produced no assembly"),
    ("LOW_COVERAGE", "reads present but below the coverage gate (or below --pathogen_min_reads)"),
    ("DETECTED_ONLY", "detect-only species with Kraken2 reads (no reconstruction)"),
    ("LOW_DEPTH", "sample below --min_reads after host removal: reconstruction skipped"),
    ("NOT_DETECTED", "no Kraken2 read in the species clade"),
]
TAXID_RE = re.compile(r"_kraken:taxid\|(\d+)")
READ_DB = {"card": "amr", "vfdb": "vf", "plasmidfinder": "mge"}


class ContractError(Exception):
    pass


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def by_suffix(directory, suffix):
    found = {}
    if os.path.isdir(directory):
        for name in sorted(os.listdir(directory)):
            if name.endswith(suffix):
                found[name[: -len(suffix)]] = os.path.join(directory, name)
    return found


def col(row, *names):
    for n in names:
        if row.get(n) is not None:
            return row[n]
    return ""


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_mlst(path):
    """(status, scheme, st, alleles) from mlst's tab line: FILE SCHEME ST locus(allele)..."""
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    if not lines:
        raise ContractError("%s is empty" % path)
    f = lines[0].split("\t")
    scheme = f[1] if len(f) > 1 else "-"
    st = f[2] if len(f) > 2 else "-"
    alleles = f[3:]
    if scheme == "-":
        return "NO_SCHEME", "", "", ""
    exact = bool(alleles) and all(re.fullmatch(r"[^()\s]+\(\d+\)", a) for a in alleles)
    return ("TYPED" if st.isdigit() and exact else "PARTIAL"), scheme, (st if st.isdigit() else ""), ";".join(alleles)


def parse_checkm2(paths):
    out = {}
    for p in paths:
        for r in read_tsv(p):
            out[r["Name"]] = (r.get("Completeness", ""), r.get("Contamination", ""))
    return out


def parse_quast(paths):
    out = {}
    for p in paths:
        rows = [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
        if not rows or rows[0][0] != "Assembly":
            raise ContractError("%s is not a QUAST report.tsv" % p)
        metrics = {r[0]: r[1:] for r in rows[1:]}
        for i, name in enumerate(rows[0][1:]):
            pick = {k: (metrics[m][i] if m in metrics and i < len(metrics[m]) else "")
                    for k, m in (("contigs", "# contigs"), ("total_length", "Total length"), ("n50", "N50"))}
            out[name] = pick
    return out


def parse_amrfinder(path):
    return sorted({col(r, "Element symbol", "Gene symbol") for r in read_tsv(path)
                   if col(r, "Type", "Element type") == "AMR" and col(r, "Element symbol", "Gene symbol")})


def parse_abricate_tab(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
    if not rows:
        return []
    if not rows[0][0].startswith("#FILE"):
        raise ContractError("%s is not an ABRicate tab report" % path)
    gi = rows[0].index("GENE")
    return sorted({r[gi] for r in rows[1:] if len(r) > gi})


def parse_mob(path):
    reps = set()
    for r in read_tsv(path):
        for rep in col(r, "rep_type(s)").split(","):
            if rep.strip() and rep.strip() != "-":
                reps.add("mob:" + rep.strip())
    return sorted(reps)


def parse_read_abricate(path):
    """{taxid: {'amr': set, 'vf': set, 'mge': set}} from a read-level ABRicate CSV (FUNCTIONAL_READS)."""
    out = defaultdict(lambda: {"amr": set(), "vf": set(), "mge": set()})
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return out
    if not rows[0] or not rows[0][0].startswith("#FILE"):
        raise ContractError("%s is not an ABRicate CSV" % path)
    si, gi, di = rows[0].index("SEQUENCE"), rows[0].index("GENE"), rows[0].index("DATABASE")
    for r in rows[1:]:
        if len(r) <= max(si, gi, di):
            continue
        m, kind = TAXID_RE.search(r[si]), READ_DB.get(r[di])
        if m and kind:
            out[m.group(1)][kind].add(r[gi])
    return out


class Results:
    def __init__(self, root):
        def d(name):
            return os.path.join(str(root), name)
        self.gate = by_suffix(d("gate"), ".gate.tsv")
        self.species = by_suffix(d("species"), ".species.tsv")
        self.assembly = by_suffix(d("assemblies"), ".fasta")
        self.mlst = by_suffix(d("mlst"), ".tsv")
        self.amr = by_suffix(d("amrfinderplus"), ".tsv")
        self.vfdb = by_suffix(d("abricate"), ".vfdb.txt")
        self.plasmid = by_suffix(d("abricate"), ".plasmidfinder.txt")
        self.mob = by_suffix(d("mobsuite"), ".mobtyper_results.txt")
        self.mixed = by_suffix(d("mixed"), ".mixed.tsv")
        self.checkm2 = parse_checkm2(by_suffix(d("checkm2"), "_checkm2_report.tsv").values())
        self.quast = parse_quast(by_suffix(d("quast"), "_quast.tsv").values())
        self.reads = {s: parse_read_abricate(p) for s, p in by_suffix(d("functional"), "_abricate.csv").items()}


def add_note(d, text):
    d["notes"] = (d["notes"] + "; " + text) if d.get("notes") else text


def evaluate_unit(unit, level, res):
    """(status, details) of one reconstruction unit that was extracted and gated."""
    d = {}
    if unit not in res.gate:
        # extracted/pooled but lost before the coverage gate (no mapping output): reported, never dropped
        add_note(d, "no coverage-gate result (reads did not map or the unit was not attempted)")
        return "LOW_COVERAGE", d
    g = read_tsv(res.gate[unit])[0]
    d["depth"], d["breadth"] = g["mean_depth"], g["breadth"]
    if g["pass"] != "true":
        add_note(d, "below the coverage gate")
        return "LOW_COVERAGE", d
    fasta = res.assembly.get(unit)
    if not fasta or os.path.getsize(fasta) == 0:
        add_note(d, "Flye/Medaka produced no assembly (see pipeline_info/execution_trace_*.txt)")
        return "ASSEMBLY_FAILED", d
    d["completeness"], d["contamination"] = res.checkm2.get(unit, ("", ""))
    d.update(res.quast.get(unit, {}))
    if unit not in res.species:
        raise ContractError("no species confirmation for assembled unit %s" % unit)
    s = read_tsv(res.species[unit])[0]
    d["ani"] = s["ani"]
    d["ani_reference"] = "%s (%s)" % (s["best_species"], s["best_accession"]) if s["best_species"] else ""
    if s["confirmed"] != "true":
        add_note(d, "species not confirmed (%s): no typing is reported under the panel species name" % s["reason"])
        return "SPECIES_NOT_CONFIRMED", d
    if level == "patient":
        if unit not in res.mixed:
            raise ContractError("no mixed-strain check for patient-level unit %s" % unit)
        m = read_tsv(res.mixed[unit])[0]
        d["mixed_het_per_kb"] = m["het_per_kb"]
        if m["verdict"] == "mixed":
            add_note(d, "%s minor-allele sites per kb (> %s): more than one strain; no ST reported"
                     % (m["het_per_kb"], m["max_het_per_kb"]))
            return "MIXED_STRAINS", d
        if m["verdict"] == "not_evaluable":
            add_note(d, "mixed-strain check not evaluable (%s sites with enough depth)" % m["evaluable_sites"])
    if unit not in res.mlst:
        raise ContractError("no mlst result for typable unit %s" % unit)
    status, d["scheme"], d["st"], d["alleles"] = parse_mlst(res.mlst[unit])
    amr = parse_amrfinder(res.amr[unit]) if unit in res.amr else None
    vf = parse_abricate_tab(res.vfdb[unit]) if unit in res.vfdb else None
    mge = (parse_abricate_tab(res.plasmid[unit]) if unit in res.plasmid else []) + \
          (parse_mob(res.mob[unit]) if unit in res.mob else [])
    d.update(amr_genes=";".join(amr or []), vf_genes=";".join(vf or []), mge_genes=";".join(mge), gene_source="contig")
    not_run = [name for name, value in (("AMRFinderPlus", amr), ("ABRicate VFDB", vf)) if value is None]
    if not_run:
        add_note(d, "not run: " + ", ".join(not_run))
    return status, d


def read_level_genes(sample, taxids, res):
    hits = res.reads.get(sample, {})
    amr, vf, mge = set(), set(), set()
    for t in taxids:
        if t in hits:
            amr |= hits[t]["amr"]
            vf |= hits[t]["vf"]
            mge |= hits[t]["mge"]
    if not (amr or vf or mge):
        return {}
    return {"amr_genes": ";".join(sorted(amr)), "vf_genes": ";".join(sorted(vf)),
            "mge_genes": ";".join(sorted(mge)), "gene_source": "read-level"}


def base_row(p, level, sample, patient, contributing, reads):
    row = {c: "" for c in COLUMNS}
    row.update(level=level, sample=sample, patient=patient, contributing_samples=contributing, species=p["species"],
               tier=p["tier"], who_category=p["who_category"], action=p["action"], taxid=p["taxid"], reads=reads)
    return row


def build_report(panel_path, cand_path, samples_path, units_path, results_root):
    panel = {p["species"]: p for p in read_tsv(panel_path)}
    cands = {(c["sample"], c["species"]): c for c in read_tsv(cand_path)}
    samples = {s["sample"]: s for s in read_tsv(samples_path)}
    units = read_tsv(units_path)
    res = Results(results_root)
    patient_rows, pooled = [], {}
    for u in units:
        if u["level"] != "patient":
            continue
        p = panel[u["species"]]
        status, d = evaluate_unit(u["unit_id"], "patient", res)
        contributing = u["samples"].split(";")
        reads = sum(int(cands[(s, u["species"])]["clade_reads"]) for s in contributing)
        row = base_row(p, "patient", "", u["patient"], ";".join(contributing), reads)
        row["status"] = status
        row.update(d)
        patient_rows.append(row)
        for s in contributing:
            pooled[(s, u["species"])] = "%s (pooled, patient %s)" % (status, u["patient"])
    rows = []
    for sid in sorted(samples):
        for p in panel.values():
            c = cands.get((sid, p["species"]))
            if c is None:
                raise ContractError("no candidate row for %s x %s" % (sid, p["species"]))
            clade = int(c["clade_reads"])
            row = base_row(p, "sample", sid, samples[sid]["patient"], sid, clade)
            row.update(bracken_reads=c["bracken_reads"], bracken_fraction=c["bracken_fraction"])
            taxids = p["taxids"].split(";")
            decision = c["decision"]
            if decision == "low_depth":
                row["status"] = "LOW_DEPTH"
                add_note(row, "sample below --min_reads: reconstruction skipped")
                if p["action"] == "reconstruct":
                    row.update(read_level_genes(sid, taxids, res))
            elif decision == "detect":
                row["status"] = "DETECTED_ONLY" if clade > 0 else "NOT_DETECTED"
            elif decision == "not_detected":
                row["status"] = "NOT_DETECTED"
            elif decision == "below_min_reads":
                row["status"] = "LOW_COVERAGE"
                add_note(row, "%d reads < --pathogen_min_reads: not extracted" % clade)
                row.update(read_level_genes(sid, taxids, res))
            elif decision == "extract":
                status, d = evaluate_unit("%s.%s" % (sid, p["slug"]), "sample", res)
                row["status"] = status
                row.update(d)
                if status in ("LOW_COVERAGE", "ASSEMBLY_FAILED"):
                    row.update(read_level_genes(sid, taxids, res))
            else:
                raise ContractError("unknown decision '%s' for %s x %s" % (decision, sid, p["species"]))
            row["pooled_result"] = pooled.get((sid, p["species"]), "")
            rows.append(row)
    keys = [(r["sample"], r["species"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ContractError("more than one status for a sample x species")
    return rows, patient_rows, samples


CSS = """body{font-family:Helvetica,Arial,sans-serif;margin:24px;color:#1f2328}
table{border-collapse:collapse;margin:8px 0 24px;font-size:13px}
th,td{border:1px solid #d0d7de;padding:4px 8px;text-align:left;vertical-align:top}th{background:#f6f8fa}
.st{font-weight:700;padding:1px 6px;border-radius:4px;white-space:nowrap}
.TYPED{background:#d1f0d9}.PARTIAL,.NO_SCHEME{background:#fff3c4}
.MIXED_STRAINS,.SPECIES_NOT_CONFIRMED,.ASSEMBLY_FAILED{background:#ffd8d3}
.LOW_COVERAGE,.LOW_DEPTH,.DETECTED_ONLY{background:#e7ecf3}.NOT_DETECTED{color:#57606a}
.note{color:#57606a;font-size:12px}"""
TABLE_COLS = [("species", "Species"), ("status", "Status"), ("reads", "Reads"), ("bracken_fraction", "Bracken fraction"),
              ("depth", "Depth"), ("breadth", "Breadth"), ("completeness", "CheckM2 compl. %"),
              ("contamination", "CheckM2 contam. %"), ("ani", "ANI %"), ("ani_reference", "Best ANI hit"),
              ("scheme", "MLST scheme"), ("st", "ST"), ("alleles", "Alleles"), ("amr_genes", "AMR"),
              ("vf_genes", "VF"), ("mge_genes", "MGE"), ("gene_source", "Gene source"),
              ("pooled_result", "Pooled result"), ("notes", "Notes")]


def table(rows, patient=False):
    e = html.escape
    cols = ([("contributing_samples", "Contributing samples")] if patient else []) + TABLE_COLS
    out = ["<table><tr>" + "".join("<th>%s</th>" % e(t) for _k, t in cols) + "</tr>"]
    for r in rows:
        cells = []
        for k, _t in cols:
            v = str(r.get(k, ""))
            cells.append("<td><span class='st %s'>%s</span></td>" % (e(v), e(v)) if k == "status"
                         else "<td>%s</td>" % e(v.replace(";", "; ")))
        out.append("<tr>" + "".join(cells) + "</tr>")
    out.append("</table>")
    return "\n".join(out)


def render(sid, smeta, srows, prows, header):
    e = html.escape
    out = ["<!doctype html><html><head><meta charset='utf-8'><title>PNEUMA pathogens: %s</title><style>%s</style>"
           "</head><body>" % (e(sid), CSS),
           "<h1>Pathogen report: %s</h1>" % e(sid),
           "<p>Patient: %s &middot; negative control: %s &middot; reads after QC: %s &middot; Kraken2 unclassified: %s%%</p>"
           % (e(smeta["patient"] or "-"), e(smeta["is_control"]), e(smeta["n_reads"]), e(smeta["unclassified_pct"])),
           "<p class='note'>PNEUMA %s &middot; panel %s (sha256 %s) &middot; coverage gate: mean depth &ge; %s, "
           "breadth &ge; %s (provisional defaults, revisited after spike-in validation)</p>"
           % (e(header["version"]), e(header["panel_name"]), e(header["panel_sha256"][:12]),
              e(header["gate_min_depth"]), e(header["gate_min_breadth"]))]
    if header["manifest"]:
        out.append("<details><summary>Databases (database_manifest.tsv)</summary><table><tr><th>name</th><th>version"
                   "</th><th>sha256</th></tr>" + "".join("<tr><td>%s</td><td>%s</td><td>%s</td></tr>"
                                                         % (e(m["name"]), e(m["version"]), e(m["sha256"][:12]))
                                                         for m in header["manifest"]) + "</table></details>")
    for title, tiers in SECTIONS:
        out.append("<h2>%s</h2>" % e(title))
        out.append(table([r for r in srows if r["tier"] in tiers]))
    if prows:
        out.append("<h2>Patient-level (pooled) results, patient %s</h2>" % e(smeta["patient"]))
        out.append("<p class='note'>Assembled from the pooled reads of several samples of this patient because no "
                   "single sample passed the coverage gate. These are not genomes from this sample alone; this "
                   "sample's own status is in the tables above.</p>")
        out.append(table(prows, patient=True))
    out.append("<h2>Status legend</h2><table>" + "".join(
        "<tr><td><span class='st %s'>%s</span></td><td>%s</td></tr>" % (k, k, e(v)) for k, v in STATUS_HELP) + "</table>")
    out.append("<p class='note'>MGE = plasmid replicons (PlasmidFinder) and, with --run_mobsuite, MOB-suite replicon "
               "types. Read-level genes are ABRicate hits on Kraken2-tagged reads of the species clade.</p></body></html>")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--panel", "--candidates", "--samples", "--units", "--results", "--panel-file"):
        ap.add_argument(a, required=True)
    ap.add_argument("--manifest")
    ap.add_argument("--pipeline-version", default="")
    ap.add_argument("--gate-min-depth", default="")
    ap.add_argument("--gate-min-breadth", default="")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args(argv)
    try:
        rows, patient_rows, samples = build_report(a.panel, a.candidates, a.samples, a.units, a.results)
    except (ContractError, KeyError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    os.makedirs(os.path.join(a.outdir, "html"), exist_ok=True)
    with open(os.path.join(a.outdir, "pathogen_cohort.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows + patient_rows)
    header = {"version": a.pipeline_version, "panel_name": os.path.basename(a.panel_file),
              "panel_sha256": sha256(a.panel_file), "gate_min_depth": a.gate_min_depth,
              "gate_min_breadth": a.gate_min_breadth,
              "manifest": read_tsv(a.manifest) if a.manifest and os.path.exists(a.manifest) else []}
    for sid, smeta in samples.items():
        srows = [r for r in rows if r["sample"] == sid]
        prows = [r for r in patient_rows if sid in r["contributing_samples"].split(";")]
        with open(os.path.join(a.outdir, "html", "%s.pathogens.html" % sid), "w") as fh:
            fh.write(render(sid, smeta, srows, prows, header))
    counts = defaultdict(int)
    for r in rows:
        counts[r["status"]] += 1
    print("pathogen report: %d samples; %s; %d patient-level rows" % (
        len(samples), ", ".join("%s %d" % (k, counts[k]) for k in STATUSES if counts[k]), len(patient_rows)),
        file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
