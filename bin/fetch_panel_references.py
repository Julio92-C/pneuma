#!/usr/bin/env python3
"""Download the PNEUMA panel reference genomes (public NCBI RefSeq assemblies).

  fetch_panel_references.py --references panel_references.tsv --outdir panel_references

For each accession of panel_references.tsv (species, accession, role): NCBI Datasets v2 genome
package (GENOME_FASTA), written as <outdir>/<accession>.fna.gz, then
<outdir>/panel_references_index.tsv (species, accession, role, file, sha256).
Three tries per accession; any failure is fatal (exit 1).
"""
import argparse
import csv
import gzip
import hashlib
import io
import os
import sys
import time
import urllib.request
import zipfile

URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/%s/download?include_annotation_type=GENOME_FASTA"
INDEX_COLUMNS = ["species", "accession", "role", "file", "sha256"]


def http_get(accession, tries=3):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(URL % accession, timeout=300) as resp:
                return resp.read()
        except OSError:
            if attempt == tries - 1:
                raise
            time.sleep(10 * (attempt + 1))
    return b""


def fasta_from_package(accession, data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        members = [m for m in z.namelist()
                   if m.startswith("ncbi_dataset/data/%s/" % accession) and m.endswith(".fna")]
        if len(members) != 1:
            raise ValueError("%s: expected one .fna in the NCBI package, found %d" % (accession, len(members)))
        return z.read(members[0])


def run(references, outdir, getter=http_get):
    os.makedirs(str(outdir), exist_ok=True)
    with open(str(references), newline="") as fh:
        refs = list(csv.DictReader(fh, delimiter="\t"))
    rows = []
    for r in refs:
        acc = r["accession"].strip()
        fasta = fasta_from_package(acc, getter(acc))
        name = "%s.fna.gz" % acc
        with open(os.path.join(str(outdir), name), "wb") as out:
            out.write(gzip.compress(fasta, mtime=0))
        digest = hashlib.sha256(open(os.path.join(str(outdir), name), "rb").read()).hexdigest()
        rows.append({"species": r["species"].strip(), "accession": acc, "role": r["role"].strip(),
                     "file": name, "sha256": digest})
    with open(os.path.join(str(outdir), "panel_references_index.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=INDEX_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--references", required=True)
    ap.add_argument("--outdir", default="panel_references")
    args = ap.parse_args(argv)
    try:
        rows = run(args.references, args.outdir)
    except (ValueError, OSError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    print("fetched %d reference genomes" % len(rows), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
