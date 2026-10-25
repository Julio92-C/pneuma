# Julio92-C/pneuma: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3.0dev - [unreleased]

Phase 2: priority-pathogen reconstruction and the pathogen report.

### `Added`

### `Changed`

### `Fixed`

## v0.2.1 - 2026-10-25

### `Fixed`

- No MERIDIAN HTML report on SGE/Singularity compute nodes (e.g. UCL Myriad): the scheduler sets `TMPDIR` to a node-local path that Singularity does not bind, and Quarto failed with `NotFound: No such file or directory (os error 2): tmpdir`. `MERIDIAN_RUN` now exports `TMPDIR`/`TMP`/`TEMP` to a `.tmp` folder in its task directory. The read-only real-MERIDIAN test (`tests/meridian_readonly.config`) sets an unbound `TMPDIR` to reproduce it.

## v0.2.0 - 2026-10-21

### `Added`

- In-pipeline database mode: `--build_databases --db_outdir <dir>` runs only the new `PREPARE_DATABASES` subworkflow (no `--input`, no analysis). It downloads Kraken2 PlusPF (`--kraken_release`), the NCBI taxdump (`--taxdump_month`), the Hostile `human-t2t-hla` index and the ABRicate `card`, `vfdb` and `plasmidfinder` databases, each in a pinned container, publishes them to `--db_outdir` (copy) and writes `database_manifest.tsv` (`name`, `version`, `path`, `sha256`; fails on an empty sha256 or date). `--skip_kraken_download`, `--skip_taxdump_download`, `--skip_hostile_download` and `--skip_abricate_download` leave a database out.
- `--db_dir <dir>` for analysis runs: fills `--taxdump`, `--host_index`, `--abricate_db`, `--database_manifest` and `--kraken2_db` from a folder built by `--build_databases` (or `bin/prepare_databases.sh`). Explicit parameters still override; the required files are validated with clear errors.
- MultiQC now summarises the whole pipeline in the order QC -> host removal -> taxonomy -> contamination -> functional: NanoStat, SeqKit stats, nanoq, Porechop (with `--run_porechop`), Kraken2 and Bracken sections, plus custom-content tables for Hostile (reads in/out, % host removed), Recentrifuge (`pneuma_controls.rcf.stat.csv`) and ABRicate (hits per sample and database). The General Statistics table gets `Reads after QC`, `% Host removed` and `% Classified` columns. Tables are built by the new `bin/make_mqc_tables.py` (process `MULTIQC_TABLES`, also published to `multiqc/custom_tables/`) with pytest and nf-test coverage.
- `bin/build_database_manifest.py`, modules `PREPARE_KRAKEN2_DB`, `PREPARE_TAXDUMP`, `PREPARE_HOSTILE_INDEX`, `PREPARE_ABRICATE_DB`, `DATABASE_MANIFEST`, with `-stub` blocks; nf-test and pytest coverage; database mode on the metro map.

### `Changed`

- `--meridian_group` accepts a comma-separated list of metadata columns: one MERIDIAN comparison (`MERIDIAN_CONFIG` + `MERIDIAN_RUN`) per column, in parallel, with study id `<study_id>_<column>`, published to `meridian/<column>/` (a single column uses the same layout). Each column is validated (exists, at least 2 distinct non-NA levels among non-control samples) by `MERIDIAN_CONFIG`, which needs only the metadata and runs at the start of the pipeline.
- The MERIDIAN image is the single parameter `--meridian_container` (default `ghcr.io/julio92-c/meridian:v1.1.0`).
- MERIDIAN v1.1.0 stage safety net: a MERIDIAN analysis stage that fails no longer aborts the run; it is listed in `meridian/<column>/Reports/stage_errors.csv` and the report, and PNEUMA logs a Nextflow warning naming the comparison and the failed stages.
- Real-container integration test for a sparse low-biomass case (`tests/data/meridian_fixture/sparse/`): ARG hits only in samples of one group and a single VF gene; both comparisons must produce a report with no failed stages.
- Tube map redrawn: branching lines from the Kraken2 interchange, a three-in/two-out junction before MERIDIAN, MultiQC as its own branch, and filled stops marking the steps reported in MultiQC.
- `MERIDIAN_RUN` works in a writable copy of `/opt/meridian` (renv library symlinked), so the Quarto report and `Rplots.pdf` no longer fail when the image is read-only or the container runs as a non-root user.
- Real-container integration test (`tests/meridian_integration.nf.test`, tag `meridian_real`, synthetic fixture from `tests/make_meridian_fixture.py`), run on demand with `nf-test test tests/meridian_integration.nf.test -c nf-test.meridian.config --tag meridian_real --profile docker,test`.

- `--input` and `--outdir` (and `--kraken2_db`, `--taxdump`, `--host_index`) are no longer `required` in `nextflow_schema.json`; nf-schema cannot express conditional requirements, so `PIPELINE_INITIALISATION` enforces them (all required for analysis unless resolvable from `--db_dir`; `--db_outdir` plus the release parameters for `--build_databases`).
- Docs: `--build_databases` is the primary way to prepare databases; `bin/prepare_databases.sh` stays as a fallback, fetched with `curl` from the release tag (the Nextflow assets path in the old docs was wrong for newer Nextflow).

### `Fixed`

- With several `--meridian_group` columns only the first comparison ran `MERIDIAN_RUN` (the others were configured but never ran): the shared inputs (ABRicate summary, Kraken2/Bracken tables, Recentrifuge CSV) were queue channels holding one item, consumed by the first `(group, config)` tuple. They are now value channels (`.first()`), so `MERIDIAN_RUN` runs once per group. A pipeline-level nf-test (`tests/main.nf.test`, real workflow wiring) asserts one successful `MERIDIAN_RUN` per group and that `meridian/` holds only `inputs/` and the per-group folders.
- No MERIDIAN HTML report under Singularity/Apptainer: Quarto renders `templates/report.qmd` in place inside `/opt/meridian`, which is read-only there, and `run_pipeline.R` still logged "done". `MERIDIAN_RUN` now builds a writable mirror `meridian_rw/` in the task directory (symlinks to every top-level entry of `/opt/meridian`, a real copy of `templates/` and `run_pipeline.R`), points `R_PROFILE_USER` at the mirror's `.Rprofile`, verifies renv still resolves `/opt/meridian/renv/library`, redirects the Quarto/deno XDG caches into the task directory, and exits non-zero when `stages.report` is true but `Reports/pipeline_report.html` is missing. The real-MERIDIAN integration test now runs the container `--read-only` (`tests/meridian_readonly.config`) and asserts the report for each group.
- Hand MERIDIAN only per-sample Recentrifuge triplets (drop EXCLUSIVE/SHARED/Details columns): the raw Recentrifuge 2.1.2 CSV made MERIDIAN's `clean_data` fail with `left_join(): Input columns in x must be unique`. New `RCF_FOR_MERIDIAN` process (`bin/rcf_for_meridian.py`) keeps column 1 plus exactly one `cnt/una/sco` triplet per samplesheet sample (controls included) and writes `meridian/inputs/pneuma.rcf.meridian.csv`; `contamination/pneuma.rcf.data.csv` stays raw.
- Recentrifuge data CSV header (Samples/Stats labels) required by MERIDIAN v1.0.0: Recentrifuge 2.1.2 leaves the first cell of row 1 and row 2 empty; `RECENTRIFUGE` now normalises them with `sed` so that MERIDIAN's `clean_data` step can parse the CSV correctly.
- MERIDIAN renv activation under Singularity: on HPC clusters with unwritable HOME, `$HOME/.cache/R/renv/sandbox` activation failed, falling back to an empty library. `MERIDIAN_RUN` now exports `RENV_CONFIG_SANDBOX_ENABLED=FALSE` and pins the sandbox directory to `$PWD/.renv_root`.
- KrakenTools `kreport2krona` outputs were published to a stray top-level `krakentools/` folder. They are now intermediate (not published); every process in `conf/modules.config` was checked, and the same fix applies to `CAT_FASTQ` (`cat/`) and `UNTAR` (`untar/`, the unpacked Kraken2 database). `PORECHOP_ABI` now publishes only its log, to `preprocess/porechop_abi/`. An nf-test asserts the published top-level folders.
- Per-sample `functional/<sample>_abricate.csv` repeated the header once per database; it now has a single header (`summary_reportClean.csv` is unchanged).
- NanoPlot wrote an identically named `NanoStats.txt` for every sample, so MultiQC saw only one sample; the files are now named `<sample>.NanoStats.txt`. nanoq runs with `-vv` so that its report is parseable by MultiQC.
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
