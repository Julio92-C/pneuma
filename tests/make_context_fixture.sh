#!/usr/bin/env bash
# Real-data fixture for tests/context/ (tag context_real): reads with ABRicate hits from one public chicken caecal
# ONT run (SRR37217670, BioProject PRJNA1406192), in the shape FUNCTIONAL_READS hands to GENOMIC_CONTEXT:
#   tests/data/context_real/S1_kraken2ID_PFP.fasta   reads >= 3 kb, ids tagged "_kraken:taxid|0"
#   tests/data/context_real/S1_abricate.csv          ABRicate card, vfdb, plasmidfinder (80/80) on those reads
# Needs curl, gzip, awk and Docker. Usage: bash tests/make_context_fixture.sh [n_reads]   (default 4000)
set -euo pipefail
N=${1:-4000}
OUT=tests/data/context_real
URL=https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR372/070/SRR37217670/SRR37217670_1.fastq.gz
mkdir -p "$OUT"
curl -sSL "$URL" | gzip -dc | awk -v n="$N" '
    NR % 4 == 1 { id = substr($2, 1, index($2 "/", "/") - 1); if (id == "") id = substr($1, 2) }
    NR % 4 == 2 && length($0) >= 3000 { print ">" id "_kraken:taxid|0"; print; if (++k >= n) exit }' \
    > "$OUT/S1_kraken2ID_PFP.fasta" || true
docker run --rm -v "$PWD/$OUT":/w -w /w quay.io/biocontainers/abricate:1.0.1--ha8f3691_1 sh -c '
    for db in card vfdb plasmidfinder; do abricate --db $db --minid 80 --mincov 80 --csv --quiet S1_kraken2ID_PFP.fasta; done' \
    > "$OUT/S1_abricate.csv"
echo "$(grep -c '>' "$OUT/S1_kraken2ID_PFP.fasta") reads, $(grep -vc '^#FILE' "$OUT/S1_abricate.csv") ABRicate hits in $OUT"
