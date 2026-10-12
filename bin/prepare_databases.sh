#!/usr/bin/env bash
# Download and pin PNEUMA Phase-1 databases; write database_manifest.tsv.
# Usage: prepare_databases.sh --outdir DIR --kraken-release YYYYMMDD --taxdump-month YYYY-MM
#                             [--hostile-image IMG] [--abricate-image IMG] [--dry-run]
# Container runtime: singularity/apptainer if available, else docker (override with PNEUMA_CONTAINER_RT).
set -euo pipefail

OUTDIR="" ; KRELEASE="" ; TMONTH="" ; DRY=0
HOSTILE_IMAGE="docker://quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0"   # same tag as modules/local/hostile_clean
ABRICATE_IMAGE="docker://quay.io/biocontainers/abricate:1.0.1--ha8f3691_1"   # same image as modules/local/abricate_reads
ABRICATE_DBS="card vfdb plasmidfinder"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)         OUTDIR="$2"; shift 2 ;;
    --kraken-release) KRELEASE="$2"; shift 2 ;;
    --taxdump-month)  TMONTH="$2"; shift 2 ;;
    --hostile-image)  HOSTILE_IMAGE="$2"; shift 2 ;;
    --abricate-image) ABRICATE_IMAGE="$2"; shift 2 ;;
    --dry-run)        DRY=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$OUTDIR"   ]] || { echo "--outdir is required" >&2; exit 2; }
[[ -n "$KRELEASE" ]] || { echo "--kraken-release YYYYMMDD is required (see https://benlangmead.github.io/aws-indexes/k2)" >&2; exit 2; }
[[ -n "$TMONTH"   ]] || { echo "--taxdump-month YYYY-MM is required (match the Kraken2 release month)" >&2; exit 2; }

K_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_${KRELEASE}.tar.gz"
T_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${TMONTH}-01.zip"

# Pick the container runtime.
RT="${PNEUMA_CONTAINER_RT:-}"
if [[ -z "$RT" ]]; then
  if command -v singularity >/dev/null 2>&1; then RT=singularity
  elif command -v apptainer >/dev/null 2>&1; then RT=apptainer
  elif command -v docker >/dev/null 2>&1; then RT=docker
  elif [[ $DRY -eq 1 ]]; then RT=singularity
  else echo "no container runtime found (need singularity, apptainer or docker)" >&2; exit 2
  fi
fi

if [[ $DRY -eq 0 ]]; then
  mkdir -p "$OUTDIR"
  OUTDIR=$(cd "$OUTDIR" && pwd)
fi

# ctr IMAGE CMD... : print the command line that runs CMD in IMAGE with OUTDIR writable
# and HOSTILE_CACHE_DIR set. Every path/argument is shell-quoted with printf %q because the
# line is later eval'd (so an OUTDIR containing spaces works).
ctr() {
  local image="$1"; shift
  local qout qhost qargs
  qout=$(printf '%q' "$OUTDIR")
  qhost=$(printf '%q' "$OUTDIR/hostile")
  qargs=$(printf ' %q' "$@")
  if [[ "$RT" == "docker" ]]; then
    echo "docker run --rm -u \$(id -u):\$(id -g) -e HOME=/tmp -e HOSTILE_CACHE_DIR=$qhost -v $qout:$qout ${image#docker://}$qargs"
  else
    echo "HOSTILE_CACHE_DIR=$qhost $RT exec --bind $qout $image$qargs"
  fi
}

run() { if [[ $DRY -eq 1 ]]; then echo "$*"; else eval "$@"; fi; }
sha() { if [[ $DRY -eq 1 ]]; then echo "dry-run"; else sha256sum "$1" | cut -d' ' -f1; fi; }

run "mkdir -p '$OUTDIR/kraken2_pluspf_${KRELEASE}' '$OUTDIR/taxdump_${TMONTH}' '$OUTDIR/hostile' '$OUTDIR/abricate'"
run "wget -nv -O '$OUTDIR/k2_pluspf_${KRELEASE}.tar.gz' '$K_URL'"
run "tar -xzf '$OUTDIR/k2_pluspf_${KRELEASE}.tar.gz' -C '$OUTDIR/kraken2_pluspf_${KRELEASE}'"
run "wget -nv -O '$OUTDIR/taxdmp_${TMONTH}-01.zip' '$T_URL'"
run "unzip -o -q '$OUTDIR/taxdmp_${TMONTH}-01.zip' nodes.dmp names.dmp -d '$OUTDIR/taxdump_${TMONTH}'"
run "$(ctr "$HOSTILE_IMAGE" hostile fetch --aligner minimap2 --name human-t2t-hla)"

# Current ABRicate databases, fetched with the same image the pipeline runs (format match).
for db in $ABRICATE_DBS; do
  run "$(ctr "$ABRICATE_IMAGE" abricate-get_db --db "$db" --dbdir "$OUTDIR/abricate" --force)"
done
# abricate-get_db only downloads the FASTA ("sequences"); --setupdb builds the BLAST indexes that
# `abricate --datadir` needs at search time (and registers the DB date shown by `abricate --list`).
run "$(ctr "$ABRICATE_IMAGE" abricate --setupdb --datadir "$OUTDIR/abricate")"

MANIFEST="$OUTDIR/database_manifest.tsv"
if [[ $DRY -eq 0 ]]; then
  HOSTILE_MMI=$(find "$OUTDIR/hostile" -name 'human-t2t-hla*.mmi' | head -1)
  # abricate --list: DATABASE SEQUENCES DBTYPE DATE (columns 1 and 4)
  ABR_LIST=$(eval "$(ctr "$ABRICATE_IMAGE" abricate --list --datadir "$OUTDIR/abricate")")
  [[ -n "$HOSTILE_MMI" ]] || { echo "ERROR: no human-t2t-hla*.mmi found under $OUTDIR/hostile (hostile fetch failed?)" >&2; exit 1; }
  # sha256 of FILE, failing loudly (a failure inside printf arguments would be silently swallowed).
  sha_or_die() {
    local v
    v=$(sha "$1") || v=""
    [[ -n "$v" ]] || { echo "ERROR: could not compute sha256 of $1" >&2; exit 1; }
    printf '%s' "$v"
  }
  SHA_K2=$(sha_or_die "$OUTDIR/k2_pluspf_${KRELEASE}.tar.gz")
  SHA_TX=$(sha_or_die "$OUTDIR/taxdmp_${TMONTH}-01.zip")
  SHA_HO=$(sha_or_die "$HOSTILE_MMI")
  ROWS=""
  for db in $ABRICATE_DBS; do
    adate=$(printf '%s\n' "$ABR_LIST" | awk -v d="$db" '$1==d {print $4}')
    [[ -n "$adate" ]] || { echo "ERROR: could not parse the date of ABRicate database '$db' from 'abricate --list'" >&2; exit 1; }
    asha=$(sha_or_die "$OUTDIR/abricate/$db/sequences")
    ROWS+=$(printf 'abricate_%s\t%s\t%s\t%s' "$db" "$adate" "$OUTDIR/abricate/$db" "$asha")$'\n'
  done
  {
    printf 'name\tversion\tpath\tsha256\n'
    printf 'kraken2_pluspf\t%s\t%s\t%s\n' "$KRELEASE" "$OUTDIR/kraken2_pluspf_${KRELEASE}" "$SHA_K2"
    printf 'ncbi_taxdump\t%s\t%s\t%s\n' "$TMONTH" "$OUTDIR/taxdump_${TMONTH}" "$SHA_TX"
    printf 'hostile_index\thuman-t2t-hla\t%s\t%s\n' "$HOSTILE_MMI" "$SHA_HO"
    printf '%s' "$ROWS"
  } > "$MANIFEST"
  echo "Wrote $MANIFEST"
  echo "--kraken2_db $OUTDIR/kraken2_pluspf_${KRELEASE}"
  echo "--taxdump $OUTDIR/taxdump_${TMONTH}"
  echo "--host_index $HOSTILE_MMI"
  echo "--abricate_db $OUTDIR/abricate"
  echo "--database_manifest $MANIFEST"
fi
