# PNEUMA — design spec

**Date:** 2026-10-08
**Owner:** Julio Cesar Ortega Cambara
**Status:** approved design, pre-implementation
**Home:** this spec moves into the new `Julio92-C/pneuma` repo as its first commit (it is drafted here because that repo does not exist yet).

## 1. Purpose

PNEUMA is an end-to-end Nextflow (DSL2) pipeline for long-read (ONT) shotgun metagenomics of respiratory samples (exhaled breath, sputum). It replaces the current chain of separately launched workflows and SGE scripts (`JC_2026-09_EBSP/Bash_scripts`: taxprofiler, wf-metagenomics, standalone Bracken, Recentrifuge, ABRicate, kractor + nf-core/bacass, merge/samplesheet helpers) with one run that:

1. processes every sample from raw FASTQ with one pinned toolset and database set;
2. feeds MERIDIAN (`ghcr.io/julio92-c/meridian`) for cohort-level ecology and resistome/virulome/mobilome statistics;
3. reconstructs, types (MLST) and annotates priority respiratory pathogens and commensal AMR/VF reservoirs, producing a per-sample pathogen report.

**v1 use case:** cohort research pipeline (batch of samples on UCL Myriad). Turnaround speed is not a v1 goal.

**Primary dataset:** the two-year lung cohort, pooled: 2025 UCLH (3 runs, Apr/May/Jun) + 2026 RFH (4 libraries EBSP0203/0304/0405R/0506R, multiple flow cells). The point of pooling through PNEUMA is to remove the phase-vs-pipeline confounding recorded in the dissertation Methods §3.10.

## 2. Non-goals (deferred)

- minimap2 classifier / custom databases (e.g. a genome-level eHOMD build). The Taxonomy subworkflow keeps a classifier interface so this slots in later. The existing `eHOMD_wf16s_db` is a 16S-only amplicon DB and is not suitable for shotgun classification.
- Whole-metagenome assembly + binning.
- TB-Profiler, pneumococcal serotyping, *H. influenzae* capsule typing, Kleborate.
- Fast per-sample clinical mode (6 h target).
- RNA viruses.
- `contamscope` as a second contamination module.
- Importing pre-existing Kraken/Bracken reports (raw FASTQ only).

## 3. Repository

- Name: **PNEUMA**; repo `Julio92-C/pneuma`, **private**, **MIT license** from day one. Visibility stays private until the UWL KE team confirms the 2025 invention disclosure (IDF) does not restrict release.
- Scaffold: nf-core pipeline template (Nextflow ≥ 24.10, nf-schema, nf-test). Not submitted to nf-core.
- MERIDIAN stays a separate repo and is consumed as a pinned container.

```
pneuma/
├── main.nf · nextflow.config · nextflow_schema.json · LICENSE (MIT)
├── workflows/pneuma.nf
├── subworkflows/local/   preprocess · taxonomy · contamination · functional_reads
│                         · pathogen_reconstruction · pathogen_report
├── modules/nf-core/      nanoq, porechop/abi, hostile, nanoplot, seqkit/stats, kraken2,
│                         bracken, krakentools, krona, abricate, minimap2, flye, medaka,
│                         checkm2, quast, skani, bakta, amrfinderplus, mlst, mobsuite, multiqc
├── modules/local/        cat_runs · tag_reads_taxid · resolve_panel · coverage_gate
│                         · recentrifuge · pathogen_report · meridian
├── bin/                  Python helpers for local modules
├── assets/               schema_input.json · pathogen_panel.csv · panel_references.tsv
├── conf/                 base · modules · myriad · test
├── bin/prepare_databases.sh  DB + reference download, writes database_manifest.tsv with sha256
├── tests/                nf-test
└── docs/ · CITATIONS.md · CHANGELOG.md
```

## 4. Inputs

### 4.1 Samplesheet (`--input`, validated by `assets/schema_input.json`)

| column | required | notes |
|---|---|---|
| `sample` | yes | canonical sample ID, regex `^[A-Za-z0-9][A-Za-z0-9-]*$` (no `_` or `.`: MERIDIAN v1.0.0 splits on them); repeated across rows when a sample was sequenced in >1 run/flow cell |
| `run_id` | yes | flow cell / run identifier |
| `batch` | yes | e.g. `UCLH_2025`, `RFH_2026`; passed to MERIDIAN as covariate |
| `site` | no | hospital |
| `fastq` | yes | one merged `.fastq.gz` per sample×run (barcode files already concatenated) |
| `is_control` | yes | `true` for negative controls |

Rows with the same `sample` are concatenated at read level before QC (replaces count-summing in `combine_batches.R` / `prep_inputs.R`). Re-extractions that the analyst wants pooled (e.g. EBR0-2466 → EBF0-2466) are expressed by giving them the same `sample`.

### 4.2 Pathogen panel (`--pathogen_panel`, default `assets/pathogen_panel.csv`)

Columns: `species, tier, who_category, action, reference_accession, confusable_with, ani_min`.
`action ∈ {reconstruct, detect}`. Species are listed by **name**; `resolve_panel` maps names → taxids against the loaded taxonomy and fails the run with closest-name suggestions if a name does not resolve.

| tier | species | WHO | action |
|---|---|---|---|
| 1 | *Streptococcus pneumoniae* | BPPL 2024 medium | reconstruct |
| 1 | *Haemophilus influenzae* | BPPL 2024 medium | reconstruct |
| 2 | *Klebsiella pneumoniae* | critical | reconstruct |
| 2 | *Acinetobacter baumannii* | critical | reconstruct |
| 2 | *Pseudomonas aeruginosa* | high | reconstruct |
| 2 | *Staphylococcus aureus* | high | reconstruct |
| 2 | *Escherichia coli* | critical | reconstruct |
| 3 | *Stenotrophomonas maltophilia* | — | reconstruct |
| 3 | *Burkholderia cenocepacia*, *Burkholderia multivorans* (Bcc; one row each) | — | reconstruct |
| 3 | *Mycobacteroides abscessus* | — | reconstruct |
| 3 | *Moraxella catarrhalis* | — | reconstruct |
| 3 | *Neisseria meningitidis* | — | reconstruct |
| D | *Mycobacterium tuberculosis* | critical | detect |
| D | *Aspergillus fumigatus*, *Candida albicans*, *Cryptococcus neoformans* | FPPL 2022 critical | detect |
| D | *Pneumocystis jirovecii* | FPPL 2022 medium | detect |
| R | *Haemophilus parainfluenzae*, *Streptococcus mitis*, *Streptococcus oralis*, *Prevotella melaninogenica*, *Veillonella parvula*, *Neisseria subflava*, *Staphylococcus epidermidis* | — (commensal AMR/VF reservoirs) | reconstruct |

`confusable_with` drives skani species confirmation (e.g. *S. pneumoniae* ↔ *S. mitis/oralis/pseudopneumoniae*; *H. influenzae* ↔ *H. haemolyticus*; *K. pneumoniae* ↔ *K. variicola/quasipneumoniae*; *S. aureus* ↔ *S. argenteus*; *E. coli* ↔ *Shigella*). Reference genomes for panel species and their confusables are listed in `assets/panel_references.tsv` and downloaded by `bin/prepare_databases.sh`.

## 5. Subworkflows and tools

### 5.1 PREPROCESS
`cat_runs` (per sample) → Porechop_ABI (optional, `--run_porechop`, default off: Dorado already trims) → Nanoq (`--min_length`, `--min_qscore`; defaults set during regression after checking the inherited taxprofiler values against short MDA-derived reads) → Hostile (minimap2 backend, `human-t2t-hla` index) → NanoPlot + seqkit stats. Samples with fewer than `--min_reads` clean reads are flagged `LOW_DEPTH`. Host reads are never published.

### 5.2 TAXONOMY
Kraken2 against the **full PlusPF** database (one pinned release, recorded with checksum), `--kraken2_confidence` exposed (default 0, revisited in regression) → Bracken (species, threshold 10, longest available read-length k-mer distribution; long-read caveat documented) → KrakenTools `combine_kreports.py` → Krona. Outputs: per-sample `.krk`, kreports, combined Kraken2 and Bracken tables.

### 5.3 CONTAMINATION
Two Recentrifuge invocations on the per-sample `.krk` files (flags `-s NORMA -y 25 -x 9606 -d -e CSV`, as in `recentrifuge_v2.sh`):
- **plain mode** (no `-c`) → `pneuma.rcf.data.csv`, fed to MERIDIAN (the `-c` mode adds CTRL/EXCLUSIVE columns that MERIDIAN's column indexing does not expect);
- **control mode** (`-c <n_controls>`, control `.krk` files passed first) → interactive HTML for contamination review only.

Count-level control subtraction is done by MERIDIAN (`metadata.controls` = all `is_control` samples). Per-run subtraction (Phase-1 style) is deferred. Runs with no control produce a warning in the run log, not an error.

### 5.4 FUNCTIONAL_READS
`tag_reads_taxid` appends `kraken:taxid|<id>` to each classified read name → ABRicate against CARD, VFDB, PlasmidFinder (min identity 80, min coverage 80) → summary table in MERIDIAN's `summary_reportClean.csv` layout. "MGE" = plasmid replicons (PlasmidFinder scope), stated in the report.

### 5.5 PATHOGEN_RECONSTRUCTION (per sample × panel species with `action = reconstruct`)
1. `resolve_panel` → taxids (+ children).
2. KrakenTools `extract_kraken_reads.py --include-children`.
3. `coverage_gate`: minimap2 of extracted reads to the panel reference → mean depth and breadth. Pass if depth ≥ `--gate_min_depth` and breadth ≥ `--gate_min_breadth` (defaults set from the spike-in validation, §8.3). Below gate → status `LOW_COVERAGE`, read-level ABRicate hits for that taxid are still reported, labelled read-level.
4. Flye → Medaka (bacterial model `r1041_e82_400bps_bacterial_methylation`).
5. CheckM2 (completeness/contamination) + QUAST.
6. skani ANI vs. the panel reference and its confusables. Best hit must be the panel species with ANI ≥ `ani_min` (default 95) → otherwise status `SPECIES_NOT_CONFIRMED` and no typing is reported under the panel species name.
7. Bakta (annotation) · AMRFinderPlus (`--organism` where supported) · ABRicate VFDB + PlasmidFinder on contigs · MOB-suite (optional, `--run_mobsuite`) · mlst (PubMLST schemes; missing loci → `PARTIAL`; no scheme → `NO_SCHEME`).

`action = detect` species skip steps 2–7 and report read counts / Bracken abundance only.

### 5.6 PATHOGEN_REPORT
Per-sample HTML + a cohort-wide TSV. Sections: Tier 1–3 pathogens, detect-only, Tier R commensal reservoirs. Exactly one status per sample × species:

`NOT_DETECTED` · `DETECTED_ONLY` (detect tier) · `LOW_DEPTH` (sample) · `LOW_COVERAGE` · `ASSEMBLY_FAILED` · `SPECIES_NOT_CONFIRMED` · `TYPED` · `PARTIAL` · `NO_SCHEME`

Each row carries reads, depth, breadth, CheckM2 completeness/contamination, ANI, ST/alleles, AMR/VF/MGE genes with source (contig vs read-level). Header records database/panel versions and the unclassified fraction.

### 5.7 MERIDIAN
Local module runs the pinned MERIDIAN container (`ghcr.io/julio92-c/meridian:v1.0.0`) with a generated `config.yaml` pointing at PNEUMA outputs (Bracken table, Kraken2 combined report, ABRicate summary, Recentrifuge CSV) plus the user's metadata file (`--meridian_metadata`, `--meridian_group`, optional `--meridian_random_effect`). No MERIDIAN patch is needed: v1.0.0 already supports `alpha_diversity.source: computed`, so `inputs.diversity_csv` is omitted. Because v1.0.0 derives the ABRicate sample ID by splitting at the first `_`, PNEUMA sample IDs must not contain underscores (enforced by the samplesheet schema). Optional: `--skip_meridian`.

### 5.8 Reporting/provenance
MultiQC over all QC; `software_versions.yml`; `database_manifest.tsv` (name, version/date, path, sha256) written by `bin/prepare_databases.sh` and copied into every run; panel file copied into results.

## 6. Infrastructure

- Development: WSL2 (Ubuntu) on Windows with Docker or conda; `-profile test` uses a tiny subsampled dataset and a toy Kraken2 DB.
- Production: `-profile myriad,singularity`. `conf/myriad.config` sets the SGE executor, converts total memory to per-core `-l mem`, sets `h_rt`, `tmpfs`, and a Scratch-resident `singularity.cacheDir`. Kraken2 full PlusPF is labelled for a ~100 GB node.
- Containers: biocontainers / nf-core module containers; MERIDIAN from `ghcr.io/julio92-c/meridian:<pinned>`.
- Databases and panel references are staged in advance with `bin/prepare_databases.sh` (Kraken2/Bracken full PlusPF, NCBI taxdump for Recentrifuge, Hostile index, CheckM2, Bakta, AMRFinderPlus, ABRicate DBs, MLST schemes, panel genomes).

## 7. Error handling

- Samplesheet: hard error on missing FASTQ, duplicate `sample+run_id`, invalid columns; warning when a `run_id` has no control.
- Panel: hard error when a species name does not resolve, with closest-name suggestions.
- `LOW_DEPTH` samples continue through taxonomy/MERIDIAN/report, skip reconstruction.
- Reconstruction failures never fail the run: the coverage gate pre-empts most; Flye/Medaka failures are captured as `ASSEMBLY_FAILED` via optional outputs (no blanket `errorStrategy 'ignore'`).
- Retry-with-more-memory on exit codes 104/134/137/139/140/143 for Kraken2, Flye, CheckM2, Bakta.
- `-resume` safe throughout.

## 8. Validation

1. **nf-test**: unit tests for every local module and `bin/` helper; pipeline-level test on `-profile test`; GitHub Actions CI.
2. **Regression vs existing results**: raw-data inventory for all 2025 + 2026 runs first (2025 Apr/May raw FASTQ are not on the laptop; locate on Myriad). Re-run and compare to `Lung_microbiome/metaOmics/Reports` (2025) and `JC_2026-09_EBSP/Reports` (2026): per-sample Bray–Curtis old vs new, ABRicate hit concordance. Differences from PlusPF-8 → full PlusPF are expected and documented, not eliminated.
3. **Spike-in ground truth** for reconstruction: public ONT isolate reads of *S. pneumoniae*, *H. influenzae*, *K. pneumoniae* spiked at 0.1–10 % into a host-depleted breath sample. Measures detection limit, sets gate defaults, checks ST accuracy. Negative test: *S. mitis*-only spike must not yield an *S. pneumoniae* ST (`SPECIES_NOT_CONFIRMED`).
4. **Full two-year cohort run** on Myriad as first production use.

## 9. Success criteria (v1)

- One command reproduces taxonomy, contamination, functional and MERIDIAN outputs for the pooled cohort from raw FASTQ, on Myriad, with `-resume`.
- No step of the legacy script chain is still needed.
- Pathogen report gives exactly one status per sample × panel species; spike-in negative test passes.
- Every run records tool versions, database versions/checksums and panel version.
