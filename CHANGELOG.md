# Julio92-C/pneuma: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
