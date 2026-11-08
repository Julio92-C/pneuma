#!/usr/bin/env python3
"""Generate the PUBLIC-DATA fixture for PNEUMA pathogen reconstruction (Phase 2).

Inputs are public NCBI RefSeq assemblies and NCBI Taxonomy plus seeded Badread simulation, so
there is no patient or real-sample data. Sample and patient IDs (PA1-S1, PD1, ...) are invented.

tests/data/pathogen_fixture/
  lineages.json        committed  NCBI lineages of panel / confusable / fixture taxa
  taxdump/             committed  tests/data/taxdump + lineages.json (deterministic)
  genomes.sha256       committed  sha256 of each genome FASTA (written on first fetch, checked after)
  expected_status.tsv  committed  truth table for tests/pathogen/pathogen_fixture.nf.test
  expected_st.tsv      committed  mlst ST of each source genome
  genomes/ reads/ samples/        generated, git-ignored (~200 Mb of simulated reads)
tests/data/pathogen_mini/         committed 120-kb TIGR4-segment fixture for the default CI suite

Subcommands (run from the repo root; `simulate` and `mini` need Docker):
  lineages  fetch lineages.json (network)          taxdump   rebuild taxdump/ (offline)
  fetch     download genomes/ (network)            simulate  Badread per DESIGN -> reads/
  kraken    reads/ -> samples/ (offline)           mini      -> tests/data/pathogen_mini/
  all       fetch + simulate + kraken
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import random
import subprocess
import sys
import time
import urllib.request
import zipfile
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIX = HERE / "data" / "pathogen_fixture"
MINI = HERE / "data" / "pathogen_mini"
BASE_TAXDUMP = HERE / "data" / "taxdump"
BADREAD = "quay.io/biocontainers/badread:0.4.2--pyhdfd78af_0"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
SEED = 20261020

# key: (public RefSeq assembly, true species taxid)
GENOMES = {
    "spn_tigr4": ("GCF_000006885.1", "1313"),     # Streptococcus pneumoniae TIGR4
    "spn_700669": ("GCF_000026665.1", "1313"),    # Streptococcus pneumoniae ATCC 700669 (Spain 23F-1)
    "hin_rd": ("GCF_000027305.1", "727"),         # Haemophilus influenzae Rd KW20
    "smi_b6": ("GCF_000027165.1", "28037"),       # Streptococcus mitis B6
    "vpa_dsm2008": ("GCF_000024945.1", "29466"),  # Veillonella parvula DSM 2008
    "rmu_dy18": ("GCF_000011025.1", "43675"),     # Rothia mucilaginosa DY-18 (non-panel background)
}
STRAIN = {"spn_tigr4": "170187"}
# How Kraken2 labels each genome's reads: species / strain / genus / unclassified ('U') or an
# explicit taxid. S. mitis reads are mostly mislabelled S. pneumoniae: the documented Kraken2
# confusion that SPECIES_NOT_CONFIRMED exists for.
LABELS = {
    "default": [("species", 0.80), ("genus", 0.15), ("U", 0.05)],
    "spn_tigr4": [("species", 0.70), ("strain", 0.10), ("genus", 0.15), ("U", 0.05)],
    "smi_b6": [("1313", 0.65), ("species", 0.25), ("genus", 0.05), ("U", 0.05)],
}
BG = [("vpa_dsm2008", 0.5), ("rmu_dy18", 0.3)]
DESIGN = [  # sample, patient, is_control, [(genome, depth)]
    ("PA1-S1", "PA1", False, [("spn_tigr4", 25)] + BG),
    ("PB1-S1", "PB1", False, [("spn_tigr4", 1)] + BG),
    ("PC1-S1", "PC1", False, [("smi_b6", 30)] + BG),
    ("PD1-S1", "PD1", False, [("hin_rd", 9)] + BG),
    ("PD1-S2", "PD1", False, [("hin_rd", 9)] + BG),
    ("PE1-S1", "PE1", False, [("spn_tigr4", 7)] + BG),
    ("PE1-S2", "PE1", False, [("spn_700669", 7)] + BG),
    ("NEG1", "", True, [("rmu_dy18", 0.1)]),
]
MINI_REF = "MINI-TIGR4-SEG"
MINI_SEGMENT = (300000, 420000)       # 120 kb of the TIGR4 chromosome
MINI_DESIGN = [("MINI-S1", "MA", 12), ("MINI-S2", "MB", 4), ("MINI-S3", "MB", 4)]

PANEL_TAXIDS = [1313, 727, 573, 470, 287, 1280, 562, 40324, 95486, 87883, 36809, 480, 487,
                1773, 746128, 5476, 5207, 42068, 729, 28037, 1303, 28132, 29466, 28449, 1282]
CONFUSABLE_TAXIDS = [257758, 726, 244366, 1463165, 985002, 623, 624]
EXTRA_TAXIDS = [43675, 170187]
RANK_ORDER = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
DMP_RANK = {"domain": "superkingdom"}   # rank name used by the Phase 1 dump / Recentrifuge 2.1.2
KRAKEN_CODE = {"superkingdom": "D", "kingdom": "K", "phylum": "P", "class": "C", "order": "O",
               "family": "F", "genus": "G", "species": "S"}
BRACKEN_HEADER = ("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
                  "new_est_reads\tfraction_total_reads")
ABRICATE_HEADER = ("#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,"
                   "%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE")


def http_get(url, tries=3):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as resp:
                return resp.read()
        except OSError:
            if attempt == tries - 1:
                raise
            time.sleep(5 * (attempt + 1))
    return b""


def gz_write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode("ascii"), mtime=0))   # mtime=0: byte-identical reruns


# ------------------------------------------------------------------ lineages and taxdump
def cmd_lineages(_args):
    ids = sorted(set(PANEL_TAXIDS + CONFUSABLE_TAXIDS + EXTRA_TAXIDS))
    out = []
    for i in range(0, len(ids), 20):
        chunk = ",".join(str(t) for t in ids[i:i + 20])
        doc = json.loads(http_get("%s/taxonomy/taxon/%s/dataset_report" % (DATASETS, chunk)))
        for rep in doc.get("reports", []):
            tax = rep["taxonomy"]
            cls = tax.get("classification", {})
            out.append({
                "taxid": str(tax["tax_id"]),
                "rank": tax.get("rank", "").lower(),
                "name": tax["current_scientific_name"]["name"],
                "basionym": tax["current_scientific_name"].get("basionym", {}).get("name", ""),
                "lineage": [[r, str(cls[r]["id"]), cls[r]["name"]] for r in RANK_ORDER if r in cls],
            })
    missing = sorted(set(str(t) for t in ids) - {o["taxid"] for o in out}, key=int)
    if missing:
        sys.exit("NCBI returned no lineage for taxid(s): " + ", ".join(missing))
    out.sort(key=lambda o: int(o["taxid"]))
    FIX.mkdir(parents=True, exist_ok=True)
    (FIX / "lineages.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print("wrote %s (%d taxa)" % (FIX / "lineages.json", len(out)))


def load_lineages():
    return json.loads((FIX / "lineages.json").read_text())


def _dmp(line):
    return [f.strip() for f in line.rstrip("\n").split("|")]


def build_taxdump(lineages, base_dir):
    """(nodes.dmp lines, names.dmp lines): the base dump plus every lineage, sorted by taxid."""
    nodes, sci, syn = {}, {}, {}
    for line in (base_dir / "nodes.dmp").read_text().splitlines():
        f = _dmp(line)
        if f[0]:
            nodes[f[0]] = (f[1], f[2])
    for line in (base_dir / "names.dmp").read_text().splitlines():
        f = _dmp(line)
        if len(f) >= 4 and f[3] == "scientific name":
            sci[f[0]] = f[1]
    nodes.setdefault("131567", ("1", "no rank"))
    sci.setdefault("131567", "cellular organisms")
    for o in lineages:
        parent = "131567"
        for rank, taxid, name in o["lineage"]:
            nodes[taxid] = (parent, DMP_RANK.get(rank, rank))
            sci[taxid] = name
            parent = taxid
        if o["taxid"] not in sci:              # strain-level taxon (TIGR4) below its species
            nodes[o["taxid"]] = (parent, "no rank")
            sci[o["taxid"]] = o["name"]
        if o["basionym"] and o["basionym"] != o["name"]:
            syn[o["taxid"]] = o["basionym"]
    node_lines = ["%s\t|\t%s\t|\t%s\t|\t\t|" % (t, p, r)
                  for t, (p, r) in sorted(nodes.items(), key=lambda kv: int(kv[0]))]
    name_lines = []
    for taxid in sorted(sci, key=int):
        name_lines.append("%s\t|\t%s\t|\t\t|\tscientific name\t|" % (taxid, sci[taxid]))
        if taxid in syn:
            name_lines.append("%s\t|\t%s\t|\t\t|\tsynonym\t|" % (taxid, syn[taxid]))
    return node_lines, name_lines


def cmd_taxdump(_args):
    nodes, names = build_taxdump(load_lineages(), BASE_TAXDUMP)
    out = FIX / "taxdump"
    out.mkdir(parents=True, exist_ok=True)
    (out / "nodes.dmp").write_text("\n".join(nodes) + "\n")
    (out / "names.dmp").write_text("\n".join(names) + "\n")
    print("wrote %s (%d nodes)" % (out, len(nodes)))


def read_taxdump(d):
    parent, rank, name = {}, {}, {}
    for line in (Path(d) / "nodes.dmp").read_text().splitlines():
        f = _dmp(line)
        parent[f[0]], rank[f[0]] = f[1], f[2]
    for line in (Path(d) / "names.dmp").read_text().splitlines():
        f = _dmp(line)
        if len(f) >= 4 and f[3] == "scientific name":
            name[f[0]] = f[1]
    return parent, rank, name


# ------------------------------------------------------------------ Kraken2-style outputs
def clades(direct, parent):
    clade = {}
    for taxid, n in direct.items():
        node = taxid
        while True:
            clade[node] = clade.get(node, 0) + n
            if node == "1":
                break
            node = parent[node]
    return clade


def kreport(assign, parent, rank, name):
    """Kraken2 report lines from {taxid: directly assigned reads, 'U': unclassified}."""
    unclassified = assign.get("U", 0)
    direct = {t: n for t, n in assign.items() if t != "U" and n}
    clade = clades(direct, parent)
    total = unclassified + clade.get("1", 0)
    children = {}
    for taxid in clade:
        if taxid != "1":
            children.setdefault(parent[taxid], []).append(taxid)

    def pct(n):
        return "%6.2f" % (100.0 * n / total if total else 0.0)

    def code(taxid, parent_code):
        if taxid == "1":
            return "R"
        mapped = KRAKEN_CODE.get(rank.get(taxid, ""))
        if mapped:
            return mapped
        base = parent_code.rstrip("0123456789")
        return "%s%d" % (base, int(parent_code[len(base):] or 0) + 1)

    lines = ["%s\t%d\t%d\tU\t0\tunclassified" % (pct(unclassified), unclassified, unclassified)]

    def walk(taxid, depth, parent_code):
        c = code(taxid, parent_code)
        lines.append("%s\t%d\t%d\t%s\t%s\t%s%s" % (pct(clade[taxid]), clade[taxid], direct.get(taxid, 0),
                                                 c, taxid, "  " * depth, name.get(taxid, taxid)))
        for child in sorted(children.get(taxid, []), key=lambda x: (-clade[x], int(x))):
            walk(child, depth + 1, c)

    if "1" in clade:
        walk("1", 0, "R")
    return lines


def label_read(key, rng, parent):
    policy = LABELS.get(key, LABELS["default"])
    x, acc, label = rng.random(), 0.0, policy[-1][0]
    for lab, frac in policy:
        acc += frac
        if x < acc:
            label = lab
            break
    species = GENOMES[key][1]
    if label == "species":
        return species
    if label == "strain":
        return STRAIN[key]
    if label == "genus":
        return parent[species]
    return label                              # 'U' or an explicit (mislabelled) taxid


def read_fastq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            yield seq, qual


def write_sample(sample, parts, outdir, parent, rank, name):
    """parts: [(genome key, reads.fastq.gz)] -> Kraken2-style files of one sample; returns n reads."""
    rng = random.Random(zlib.crc32(sample.encode()) ^ SEED)
    fq, assignments, counts, n = [], [], {}, 0
    for key, reads in parts:
        for seq, qual in read_fastq(reads):
            n += 1
            rid = "%s-r%06d" % (sample, n)
            label = label_read(key, rng, parent)
            counts[label] = counts.get(label, 0) + 1
            kmers = max(1, len(seq) - 34)
            if label == "U":
                assignments.append("U\t%s\t0\t%d\t0:%d" % (rid, len(seq), kmers))
            else:
                assignments.append("C\t%s\t%s\t%d\t%s:%d" % (rid, label, len(seq), label, kmers))
                fq.append("@%s kraken:taxid|%s\n%s\n+\n%s\n" % (rid, label, seq, qual))
    outdir.mkdir(parents=True, exist_ok=True)
    gz_write(outdir / ("%s.classified.fastq.gz" % sample), "".join(fq))
    (outdir / ("%s.kraken2.classifiedreads.txt" % sample)).write_text("\n".join(assignments) + "\n")
    (outdir / ("%s.kraken2.report.txt" % sample)).write_text("\n".join(kreport(counts, parent, rank, name)) + "\n")
    clade = clades({t: c for t, c in counts.items() if t != "U"}, parent)
    species = sorted((t for t in clade if rank.get(t) == "species"), key=lambda t: (-clade[t], int(t)))
    total = sum(clade[t] for t in species) or 1
    rows = [BRACKEN_HEADER] + ["%s\t%s\tS\t%d\t0\t%d\t%.5f" % (name[t], t, clade[t], clade[t], clade[t] / total)
                               for t in species]
    (outdir / "bracken").mkdir(parents=True, exist_ok=True)
    (outdir / "bracken" / ("%s.tsv" % sample)).write_text("\n".join(rows) + "\n")
    (outdir / "abricate").mkdir(parents=True, exist_ok=True)
    (outdir / "abricate" / ("%s_abricate.csv" % sample)).write_text(ABRICATE_HEADER + "\n")
    return n


def write_samples_tsv(path, rows):
    lines = ["sample\tpatient\tis_control\tlow_depth\tn_reads"]
    lines += ["%s\t%s\t%s\tfalse\t%d" % (s, p, "true" if c else "false", n) for s, p, c, n in rows]
    path.write_text("\n".join(lines) + "\n")


# ------------------------------------------------------------------ genomes and simulation
def fetch_genome(accession):
    data = http_get("%s/genome/accession/%s/download?include_annotation_type=GENOME_FASTA" % (DATASETS, accession))
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        members = [m for m in z.namelist()
                   if m.startswith("ncbi_dataset/data/%s/" % accession) and m.endswith(".fna")]
        if len(members) != 1:
            sys.exit("%s: expected one .fna in the NCBI package, found %s" % (accession, members))
        return z.read(members[0])


def cmd_fetch(_args):
    gdir = FIX / "genomes"
    gdir.mkdir(parents=True, exist_ok=True)
    shafile = FIX / "genomes.sha256"
    known = {}
    if shafile.exists():
        for line in shafile.read_text().splitlines():
            digest, fname = line.split("  ")
            known[fname] = digest
    lines = []
    for key in sorted(GENOMES):
        path = gdir / ("%s.fna" % key)
        if not path.exists():
            path.write_bytes(fetch_genome(GENOMES[key][0]))
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if path.name in known and known[path.name] != digest:
            sys.exit("%s changed since genomes.sha256 was written (NCBI re-issued it?)" % path.name)
        lines.append("%s  %s" % (digest, path.name))
    shafile.write_text("\n".join(lines) + "\n")
    print("genomes ok: %d" % len(lines))


def badread(ref_dir, ref_name, depth, seed, out_path):
    if out_path.exists():
        return                                  # cached: simulation is the slow step
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(".tmp")
    cmd = ["docker", "run", "--rm", "-u", "%d:%d" % (os.getuid(), os.getgid()), "-v", "%s:/g:ro" % ref_dir,
           BADREAD, "badread", "simulate", "--reference", "/g/%s" % ref_name, "--quantity", "%sx" % depth,
           "--seed", str(seed), "--length", "6000,4000", "--identity", "96,99,2",
           "--error_model", "nanopore2023", "--qscore_model", "nanopore2023"]
    with open(tmp, "wb") as fh:
        subprocess.run(cmd, stdout=fh, check=True)
    out_path.write_bytes(gzip.compress(tmp.read_bytes(), mtime=0))
    tmp.unlink()


def seed_for(*parts):
    return zlib.crc32(":".join(parts).encode()) % (2 ** 31)


def cmd_simulate(_args):
    for sample, _patient, _control, components in DESIGN:
        for key, depth in components:
            badread(FIX / "genomes", "%s.fna" % key, depth, seed_for(sample, key),
                    FIX / "reads" / ("%s.%s.fastq.gz" % (sample, key)))
            print("simulated %s %s %sx" % (sample, key, depth))


def cmd_kraken(_args):
    parent, rank, name = read_taxdump(FIX / "taxdump")
    out = FIX / "samples"
    rows = []
    for sample, patient, control, components in DESIGN:
        parts = [(key, FIX / "reads" / ("%s.%s.fastq.gz" % (sample, key))) for key, _depth in components]
        rows.append((sample, patient, control, write_sample(sample, parts, out, parent, rank, name)))
    write_samples_tsv(out / "samples.tsv", rows)
    print("wrote %s" % out)


def cmd_mini(_args):
    parent, rank, name = read_taxdump(FIX / "taxdump")
    record = (FIX / "genomes" / "spn_tigr4.fna").read_text().split(">")[1]
    seq = "".join(record.splitlines()[1:])[MINI_SEGMENT[0]:MINI_SEGMENT[1]]
    fasta = ">%s TIGR4 GCF_000006885.1 chromosome %d-%d\n" % (MINI_REF, MINI_SEGMENT[0] + 1, MINI_SEGMENT[1])
    fasta += "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n"
    seg_dir = FIX / "reads" / "mini_ref"
    seg_dir.mkdir(parents=True, exist_ok=True)
    (seg_dir / ("%s.fna" % MINI_REF)).write_text(fasta)
    gz_write(MINI / "ref" / ("%s.fna.gz" % MINI_REF), fasta)
    rows = []
    for sample, patient, depth in MINI_DESIGN:
        reads = FIX / "reads" / ("%s.mini.fastq.gz" % sample)
        badread(seg_dir, "%s.fna" % MINI_REF, depth, seed_for(sample, "mini"), reads)
        rows.append((sample, patient, False, write_sample(sample, [("spn_tigr4", reads)], MINI, parent, rank, name)))
    write_samples_tsv(MINI / "samples.tsv", rows)
    (MINI / "panel.csv").write_text(
        "species,tier,who_category,action,reference_accession,confusable_with,ani_min,amrfinder_organism\n"
        "Streptococcus pneumoniae,1,BPPL 2024 medium,reconstruct,%s,,95,Streptococcus_pneumoniae\n" % MINI_REF)
    (MINI / "panel_references.tsv").write_text(
        "species\taccession\trole\nStreptococcus pneumoniae\t%s\tpanel\n" % MINI_REF)
    ref = MINI / "ref" / ("%s.fna.gz" % MINI_REF)
    (MINI / "ref" / "panel_references_index.tsv").write_text(
        "species\taccession\trole\tfile\tsha256\nStreptococcus pneumoniae\t%s\tpanel\t%s\t%s\n"
        % (MINI_REF, ref.name, hashlib.sha256(ref.read_bytes()).hexdigest()))
    print("wrote %s" % MINI)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["lineages", "taxdump", "fetch", "simulate", "kraken", "mini", "all"])
    args = ap.parse_args(argv)
    steps = {"lineages": [cmd_lineages], "taxdump": [cmd_taxdump], "fetch": [cmd_fetch],
             "simulate": [cmd_simulate], "kraken": [cmd_kraken], "mini": [cmd_mini],
             "all": [cmd_fetch, cmd_simulate, cmd_kraken]}[args.command]
    for step in steps:
        step(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
