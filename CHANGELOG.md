# Julio92-C/pneuma: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0 - 2026-10-13

### `Added`

- In-pipeline database mode: `--build_databases --db_outdir <dir>` runs only the new `PREPARE_DATABASES` subworkflow (no `--input`, no analysis). It downloads Kraken2 PlusPF (`--kraken_release`), the NCBI taxdump (`--taxdump_month`), the Hostile `human-t2t-hla` index and the ABRicate `card`, `vfdb` and `plasmidfinder` databases, each in a pinned container, publishes them to `--db_outdir` (copy) and writes `database_manifest.tsv` (`name`, `version`, `path`, `sha256`; fails on an empty sha256 or date). `--skip_kraken_download`, `--skip_taxdump_download`, `--skip_hostile_download` and `--skip_abricate_download` leave a database out.
- `--db_dir <dir>` for analysis runs: fills `--taxdump`, `--host_index`, `--abricate_db`, `--database_manifest` and `--kraken2_db` from a folder built by `--build_databases` (or `bin/prepare_databases.sh`). Explicit parameters still override; the required files are validated with clear errors.
- `bin/build_database_manifest.py`, modules `PREPARE_KRAKEN2_DB`, `PREPARE_TAXDUMP`, `PREPARE_HOSTILE_INDEX`, `PREPARE_ABRICATE_DB`, `DATABASE_MANIFEST`, with `-stub` blocks; nf-test and pytest coverage; database mode on the metro map.

### `Changed`

- `--input` and `--outdir` (and `--kraken2_db`, `--taxdump`, `--host_index`) are no longer `required` in `nextflow_schema.json`; nf-schema cannot express conditional requirements, so `PIPELINE_INITIALISATION` enforces them (all required for analysis unless resolvable from `--db_dir`; `--db_outdir` plus the release parameters for `--build_databases`).
- Docs: `--build_databases` is the primary way to prepare databases; `bin/prepare_databases.sh` stays as a fallback, fetched with `curl` from the release tag (the Nextflow assets path in the old docs was wrong for newer Nextflow).

### `Fixed`

- `hostile fetch` no longer exists in Hostile 2.0.2 (the pinned version); `bin/prepare_databases.sh` now calls `hostile index fetch --name human-t2t-hla --minimap2`.

### `Known issues`

- Nextflow 26.04 (strict syntax by default) passes command-line values as strings, so boolean and numeric flags fail schema validation; use `-params-file` or `NXF_SYNTAX_PARSER=v1` (documented in README and `docs/usage.md`).
- ABRicate's `vfdb` download (`www.mgc.ac.cn`) was returning HTTP 500 / timing out on 2026-10-13 when tested; `--build_databases` fails at the ABRicate step until the host recovers (the same applies to `prepare_databases.sh`).

## v0.1.1 - 2026-10-13

### `Added`

- MultiQC station on the pipeline metro map.
- README "Help & version" section, Myriad profile and Singularity cache guidance, description of the `test` profile data.

### `Fixed`

- `Unrecognized config option 'validation.*'` warnings on Nextflow 26.04: bumped `nf-schema` from 2.5.1 to 2.6.1.
- `RECENTRIFUGE` version command used an invalid `\.` escape in a Groovy string, which made the module fail to compile on Nextflow 25.10.4 (the declared minimum).

### `Dependencies`

- `nf-schema` 2.5.1 -> 2.6.1.

## v0.1.0 - 2026-10-12

v0.1.0 — Phase 1 core: preprocess, taxonomy, contamination, functional reads, MERIDIAN hand-off

Initial release of Julio92-C/pneuma, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Preprocess: run pooling, optional Porechop_ABI, nanoq, NanoPlot, seqkit stats, Hostile host removal, `LOW_DEPTH` flag.
- Taxonomy: Kraken2, Bracken, KrakenTools report combination, Krona.
- Contamination: Recentrifuge (plain and control mode).
- Functional reads: ABRicate (`card`, `vfdb`, `plasmidfinder`) on Kraken-tagged reads.
- MERIDIAN hand-off (config generation and containerised run).
- `docker`, `singularity`, `test` and `myriad` profiles; `bin/prepare_databases.sh`.

### `Fixed`

### `Dependencies`

### `Deprecated`
