#!/usr/bin/env python3
"""Build MultiQC custom-content tables for the PNEUMA stages MultiQC has no module for.

Inputs (all optional, any number of files each):
  --hostile   <sample>.hostile.json        Hostile 2.x run logs
  --kraken    <sample>.kraken2.report.txt  Kraken2 reports
  --abricate  <sample>_abricate.csv        per-sample ABRicate --csv outputs
  --rcf-stat  pneuma_controls.rcf.stat.csv Recentrifuge statistics (CSV)
  --pathogen-cohort pathogen_cohort.tsv  PNEUMA pathogen report (cohort TSV)

Outputs (written to --outdir, only when the matching input is present):
  pneuma_general_stats_mqc.tsv  columns for the General Statistics table:
                                reads after QC, % host removed, % classified
  hostile_mqc.tsv               "Host removal" section (reads in/out/removed)
  abricate_mqc.tsv              "ABRicate hits" section (hits per database)
  recentrifuge_mqc.tsv          "Recentrifuge" section (per-sample statistics)
  pathogen_status_mqc.tsv       "Pathogen reconstruction" section (status counts per sample)

Sample names are taken from the file names so that they match the names
MultiQC derives for its own modules.
"""
import argparse
import csv
import json
import os
import sys
from collections import Counter, OrderedDict

ABRICATE_DBS = ("card", "vfdb", "plasmidfinder")


def _strip(name, suffixes):
    """Drop any directory prefix (POSIX or Windows) and the first matching suffix."""
    name = name.replace("\\", "/").rsplit("/", 1)[-1]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


RCF_SUFFIXES = (".kraken2.classifiedreads.txt", ".classifiedreads.txt",
                ".kraken2.classifiedreads", ".classifiedreads")


def _warn(path, reason):
    print("WARNING: %s not parsed: %s" % (path, reason), file=sys.stderr)


def _is_stub(path):
    """True for the empty `{}` placeholder JSON written by stub runs."""
    try:
        with open(str(path)) as fh:
            return json.load(fh) == {}
    except (OSError, ValueError):
        return False


def sample_name(path, suffixes):
    return _strip(os.path.basename(str(path)), suffixes)


def _fmt(value):
    if isinstance(value, float):
        return "{:.2f}".format(value)
    return str(value)


def parse_hostile(paths):
    """{sample: {reads_in, reads_out, reads_removed, pct_host_removed}}"""
    data = OrderedDict()
    for path in sorted(str(p) for p in paths):
        with open(path) as fh:
            try:
                doc = json.load(fh)
            except ValueError as exc:
                _warn(path, "invalid JSON (%s)" % exc)
                continue
        records = [r for r in (doc if isinstance(doc, list) else [doc]) if isinstance(r, dict)]
        if not any("reads_in" in r for r in records):
            if doc != {}:
                _warn(path, "no reads_in field")
            continue
        reads_in = sum(int(r.get("reads_in", 0)) for r in records)
        reads_out = sum(int(r.get("reads_out", 0)) for r in records)
        removed = reads_in - reads_out
        pct = 100.0 * removed / reads_in if reads_in else 0.0
        data[sample_name(path, (".hostile.json", ".json"))] = OrderedDict(
            reads_in=reads_in, reads_out=reads_out, reads_removed=removed, pct_host_removed=pct)
    return data


def parse_kraken(paths):
    """{sample: {unclassified, classified, pct_classified}} from Kraken2 reports."""
    data = OrderedDict()
    for path in sorted(str(p) for p in paths):
        unclassified = classified = 0
        seen = False
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 6 or not f[1].strip().isdigit():
                    continue
                seen = True
                if f[3].strip() == "U":
                    unclassified = int(f[1])
                elif f[3].strip() == "R" and f[5].strip() == "root":
                    classified = int(f[1])
        if not seen:
            _warn(path, "empty or unrecognised Kraken2 report")
            continue
        total = unclassified + classified
        data[sample_name(path, (".kraken2.report.txt", ".report.txt", ".txt"))] = OrderedDict(
            unclassified=unclassified, classified=classified,
            pct_classified=100.0 * classified / total if total else 0.0)
    return data


def parse_abricate(paths, dbs=ABRICATE_DBS):
    """{sample: {<db>: hits..., total_hits, distinct_genes}}; empty files give zeros."""
    data = OrderedDict()
    for path in sorted(str(p) for p in paths):
        per_db = Counter()
        genes = set()
        if os.path.getsize(path) == 0:
            _warn(path, "empty file (no ABRicate header)")
            continue
        with open(path, newline="") as fh:
            header = None
            for row in csv.reader(fh):
                if not row:
                    continue
                if row[0].startswith("#FILE"):
                    header = row
                    continue
                if header is None or len(row) < len(header):
                    continue
                rec = dict(zip(header, row))
                per_db[rec.get("DATABASE", "")] += 1
                genes.add((rec.get("DATABASE", ""), rec.get("GENE", "")))
        row = OrderedDict((db, per_db.get(db, 0)) for db in dbs)
        for extra in sorted(set(per_db) - set(dbs)):
            row[extra] = per_db[extra]
        row["total_hits"] = sum(per_db.values())
        row["distinct_genes"] = len(genes)
        data[sample_name(path, ("_abricate.csv", ".csv"))] = row
    return data


def parse_rcf_stat(path):
    """{sample: {metric: value}} from the Recentrifuge stat CSV (metrics are rows)."""
    data = OrderedDict()
    with open(str(path), newline="") as fh:
        rows = [r for r in csv.reader(fh) if r]
    if len(rows) < 2:
        _warn(path, "no statistics rows")
        return data
    names = [sample_name(c, RCF_SUFFIXES) for c in rows[0][1:]]
    for name in names:
        data[name] = OrderedDict()
    for row in rows[1:]:
        # short/ragged rows (e.g. "Seqs. filtered,7344.0") only fill the leading samples
        for name, value in zip(names, row[1:]):
            try:
                data[name][row[0]] = float(value)
            except ValueError:
                continue
    return data


def _write(path, header_lines, columns, data, samples=None):
    keys = list(columns)
    with open(path, "w") as out:
        for line in header_lines:
            out.write("# " + line + "\n")
        out.write("\t".join(["Sample"] + keys) + "\n")
        for sample in (samples if samples is not None else data):
            row = data.get(sample, {})
            out.write("\t".join([sample] + [_fmt(row.get(k, "")) for k in keys]) + "\n")


def write_general_stats(path, hostile, kraken):
    samples = sorted(set(hostile) | set(kraken))
    merged = OrderedDict()
    for s in samples:
        row = OrderedDict()
        if s in hostile:
            row["reads_after_qc"] = hostile[s]["reads_out"]
            row["pct_host_removed"] = hostile[s]["pct_host_removed"]
        if s in kraken:
            row["pct_classified"] = kraken[s]["pct_classified"]
        merged[s] = row
    header = [
        "plot_type: 'generalstats'",
        "pconfig:",
        "    - reads_after_qc:",
        "        title: 'Reads after QC'",
        "        description: 'Reads left after nanoq filtering and Hostile host removal (input to taxonomy)'",
        "        scale: 'Greens'",
        "        format: '{:,.0f}'",
        "    - pct_host_removed:",
        "        title: '% Host removed'",
        "        description: 'Percentage of reads removed as host by Hostile'",
        "        suffix: '%'",
        "        min: 0",
        "        max: 100",
        "        scale: 'Oranges'",
        "        format: '{:,.2f}'",
        "    - pct_classified:",
        "        title: '% Classified'",
        "        description: 'Percentage of reads classified by Kraken2'",
        "        suffix: '%'",
        "        min: 0",
        "        max: 100",
        "        scale: 'Blues'",
        "        format: '{:,.2f}'",
    ]
    _write(path, header, ["reads_after_qc", "pct_host_removed", "pct_classified"], merged, samples)


def write_hostile(path, hostile):
    header = [
        "id: 'hostile'",
        "section_name: 'Hostile host removal'",
        "description: 'Reads entering and leaving Hostile (human-t2t-hla index by default). Reads in are those that passed nanoq filtering.'",
        "plot_type: 'table'",
        "pconfig:",
        "    id: 'hostile_table'",
        "    title: 'Hostile host removal'",
        "headers:",
        "    reads_in:",
        "        title: 'Reads in'",
        "        format: '{:,.0f}'",
        "    reads_out:",
        "        title: 'Reads out'",
        "        format: '{:,.0f}'",
        "    reads_removed:",
        "        title: 'Host reads'",
        "        format: '{:,.0f}'",
        "    pct_host_removed:",
        "        title: '% Host removed'",
        "        suffix: '%'",
        "        min: 0",
        "        max: 100",
        "        scale: 'Oranges'",
        "        format: '{:,.2f}'",
    ]
    _write(path, header, ["reads_in", "reads_out", "reads_removed", "pct_host_removed"], hostile)


def write_abricate(path, abricate):
    columns = list(ABRICATE_DBS)
    for row in abricate.values():
        for key in row:
            if key not in columns and key not in ("total_hits", "distinct_genes"):
                columns.append(key)
    columns += ["total_hits", "distinct_genes"]
    header = [
        "id: 'abricate_hits'",
        "section_name: 'ABRicate hits'",
        "description: 'ABRicate hits on Kraken-tagged reads, per sample and database (card = resistance genes, vfdb = virulence factors, plasmidfinder = plasmid replicons). Identity and coverage thresholds are set by --abricate_minid and --abricate_mincov.'",
        "plot_type: 'table'",
        "pconfig:",
        "    id: 'abricate_hits_table'",
        "    title: 'ABRicate hits per database'",
        "headers:",
    ]
    for key in columns:
        title = {"total_hits": "Total hits", "distinct_genes": "Distinct genes"}.get(key, key)
        header += ["    %s:" % key, "        title: '%s'" % title, "        format: '{:,.0f}'"]
    _write(path, header, columns, abricate)


def write_rcf(path, rcf):
    columns = []
    for row in rcf.values():
        for key in row:
            if key not in columns:
                columns.append(key)
    header = [
        "id: 'recentrifuge_stats'",
        "section_name: 'Recentrifuge'",
        "description: 'Per-sample statistics from the Recentrifuge control-mode run (pneuma_controls.rcf.stat.csv).'",
        "plot_type: 'table'",
        "pconfig:",
        "    id: 'recentrifuge_table'",
        "    title: 'Recentrifuge statistics'",
        "headers:",
    ]
    for key in columns:
        header += ["    '%s':" % key, "        title: '%s'" % key, "        format: '{:,.1f}'"]
    _write(path, header, columns, rcf)


PATHOGEN_STATUSES = ["TYPED", "PARTIAL", "NO_SCHEME", "MIXED_STRAINS", "SPECIES_NOT_CONFIRMED",
                     "ASSEMBLY_FAILED", "LOW_COVERAGE", "DETECTED_ONLY", "LOW_DEPTH"]


def parse_pathogen_cohort(path):
    """{sample: {<status>: n, typed: 'Species STn; ...', pooled: '...'}} from level=sample rows."""
    data = OrderedDict()
    with open(str(path), newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if rows and "level" not in rows[0]:
        _warn(path, "not a pathogen_cohort.tsv (no 'level' column)")
        return data
    for r in rows:
        if r["level"] != "sample":
            continue
        row = data.setdefault(r["sample"], OrderedDict([(s, 0) for s in PATHOGEN_STATUSES] + [("typed", ""), ("pooled", "")]))
        if r["status"] in row:
            row[r["status"]] += 1
        if r["status"] == "TYPED":
            row["typed"] = "; ".join(filter(None, [row["typed"], "%s ST%s" % (r["species"], r["st"])]))
        if r.get("pooled_result"):
            row["pooled"] = "; ".join(filter(None, [row["pooled"], "%s: %s" % (r["species"], r["pooled_result"])]))
    return data


def write_pathogen(path, data):
    header = [
        "id: 'pathogen_status'",
        "section_name: 'Pathogen reconstruction'",
        "description: 'Per sample: number of panel species in each status (one status per sample x species; NOT_DETECTED is not counted), typed species with their ST, and patient-level pooled results involving the sample. Per-sample HTML reports are in pathogens/report/html/.'",
        "plot_type: 'table'",
        "pconfig:",
        "    id: 'pathogen_status_table'",
        "    title: 'Pathogen status per sample'",
        "headers:",
    ]
    for key in PATHOGEN_STATUSES:
        header += ["    %s:" % key, "        title: '%s'" % key, "        format: '{:,.0f}'"]
    header += ["    typed:", "        title: 'Typed (ST)'", "    pooled:", "        title: 'Pooled (patient level)'"]
    _write(path, header, PATHOGEN_STATUSES + ["typed", "pooled"], data)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hostile", nargs="*", default=[])
    parser.add_argument("--kraken", nargs="*", default=[])
    parser.add_argument("--abricate", nargs="*", default=[])
    parser.add_argument("--rcf-stat", nargs="*", default=[])
    parser.add_argument("--pathogen-cohort", nargs="*", default=[])
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    written = []

    def out(name):
        return os.path.join(args.outdir, name)

    failed = []

    hostile = parse_hostile(args.hostile)
    if args.hostile and not hostile and not all(_is_stub(p) for p in args.hostile):
        failed.append("--hostile")
    kraken = parse_kraken(args.kraken)
    if args.kraken and not kraken:
        failed.append("--kraken")
    if hostile or kraken:
        write_general_stats(out("pneuma_general_stats_mqc.tsv"), hostile, kraken)
        written.append("pneuma_general_stats_mqc.tsv")
    if hostile:
        write_hostile(out("hostile_mqc.tsv"), hostile)
        written.append("hostile_mqc.tsv")
    if args.abricate:
        abricate = parse_abricate(args.abricate)
        if abricate:
            write_abricate(out("abricate_mqc.tsv"), abricate)
            written.append("abricate_mqc.tsv")
        else:
            failed.append("--abricate")
    rcf = parse_rcf_stat(args.rcf_stat[0]) if args.rcf_stat else {}
    if rcf:
        write_rcf(out("recentrifuge_mqc.tsv"), rcf)
        written.append("recentrifuge_mqc.tsv")
    elif args.rcf_stat:
        failed.append("--rcf-stat")
    pathogen = parse_pathogen_cohort(args.pathogen_cohort[0]) if args.pathogen_cohort else {}
    if pathogen:
        write_pathogen(out("pathogen_status_mqc.tsv"), pathogen)
        written.append("pathogen_status_mqc.tsv")
    elif args.pathogen_cohort:
        failed.append("--pathogen-cohort")
    print("wrote: " + (", ".join(written) or "nothing"), file=sys.stderr)
    if failed:
        print("ERROR: no rows could be built from: " + ", ".join(failed), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
