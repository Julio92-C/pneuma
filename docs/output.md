# Julio92-C/pneuma: Output

## Introduction

This document describes the output produced by PNEUMA. All paths are relative to the results directory given with `--outdir`. Directory names are fixed by `publishDir` rules in [`conf/modules.config`](../conf/modules.config); the rule applying to every other process is `<outdir>/<first word of the process name, lower-case>`.

Per-sample file names start with the sample name from the samplesheet (`SAMPLE-01` below is a placeholder). Reads from several runs of the same sample are pooled before analysis, so there is one set of per-sample outputs per `sample`.

## Pipeline overview

- [`preprocess/`](#preprocess): read QC, host removal
- [`taxonomy/`](#taxonomy): Kraken2, Bracken, combined reports, Krona
- [`contamination/`](#contamination): Recentrifuge
- [`functional/`](#functional): ABRicate on taxid-tagged reads
- [`meridian/`](#meridian): MERIDIAN config and results
- [`multiqc/`](#multiqc): one report for all stages
- [`pipeline_info/`](#pipeline-information)

### Preprocess

<details markdown="1">
<summary>Output files</summary>

- `preprocess/nanoq/`
  - `<sample>.nanoq.stats`: read-length and quality summary from nanoq (published `*.stats` only; also parsed by MultiQC).
- `preprocess/clean_reads/`
  - `<sample>.hostile.fastq.gz`: reads after host removal with Hostile; these feed the taxonomy stage.
  - `<sample>.hostile.json`: Hostile run log (JSON).
- `preprocess/nanoplot/<sample>/`: NanoPlot HTML report, plots and summary text files (`<sample>.NanoStats.txt` is parsed by MultiQC).
- `preprocess/seqkit/`: `seqkit stats -T -a` tables.
- `preprocess/porechop_abi/<sample>.porechop_abi.log`: adapter-trimming log, only with `--run_porechop`.

</details>

Reads are optionally adapter-trimmed (`porechop_abi`, `--run_porechop`), filtered with nanoq (`--nanoq_min_length`, `--nanoq_min_qscore`) and depleted of host reads with Hostile. A sample with fewer than `--min_reads` reads after host removal is flagged `LOW_DEPTH` in the log and skipped by Bracken. Pooled reads (`cat`), adapter-trimmed reads and the unpacked Kraken2 database are intermediate and not published.

### Taxonomy

<details markdown="1">
<summary>Output files</summary>

- `taxonomy/kraken2/`
  - `<sample>.report.txt`: Kraken2 report.
  - `<sample>.classifiedreads.txt`: per-read Kraken2 assignments (also used by Recentrifuge).
- `taxonomy/bracken/`: Bracken abundance tables and Bracken-style reports per sample (`-l`, `-t`, `-r` from `--bracken_level`, `--bracken_threshold`, `--bracken_read_length`).
- `taxonomy/`
  - `kraken2_combined_reports.txt`: KrakenTools-combined Kraken2 reports for all samples, arranged for MERIDIAN.
  - `bracken_arranged.csv`: combined Bracken table arranged for MERIDIAN.
- `taxonomy/krona/`: interactive Krona HTML charts. The KrakenTools `kreport2krona` text files that feed Krona are intermediate and not published.

</details>

Classification uses the Kraken2 PlusPF index (`--kraken2_db`, `--kraken2_confidence`); the NCBI taxdump (`--taxdump`) supplies lineages.

### Contamination

<details markdown="1">
<summary>Output files</summary>

- `contamination/`
  - `pneuma.rcf.data.csv`: Recentrifuge per-taxon table across samples.
  - `pneuma.rcf.html`: interactive Recentrifuge report (plain mode).
  - `pneuma_controls.rcf.*`: control-mode report, written only when the samplesheet marks negative controls (`is_control = true`).

</details>

### Functional

<details markdown="1">
<summary>Output files</summary>

- `functional/`
  - `<sample>_abricate.csv`: ABRicate hits for the Kraken-tagged reads of one sample, searched against the databases in `--abricate_dbs` (default `card,vfdb,plasmidfinder`) at `--abricate_minid` / `--abricate_mincov` (default 80 / 80). One header line, then one row per hit; the `DATABASE` column tells the databases apart. A sample without hits has the header only.
  - `summary_reportClean.csv`: hits of all samples merged into one table.

</details>

The tagged read FASTA files that feed ABRicate are intermediate and not published.

### MERIDIAN

<details markdown="1">
<summary>Output files</summary>

- `meridian/inputs/pneuma.rcf.meridian.csv`: the Recentrifuge table trimmed for MERIDIAN (shared by all comparisons).
- `meridian/<group>/`: one directory per column listed in `--meridian_group` (a single column gives a single directory)
  - `config.yaml`: the MERIDIAN configuration generated from `assets/meridian_config.template.yaml` (study id `<--meridian_study_id>_<group>`).
  - `Datasets/`: tables produced by MERIDIAN.
  - `Figures/`: figures produced by MERIDIAN.
  - `Reports/`: the MERIDIAN HTML report.

</details>

Only produced when `--skip_meridian` is not set. See the [MERIDIAN repository](https://github.com/Julio92-C/meridian) for how to read these files.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.
  - `custom_tables/`: the `*_mqc.tsv` custom-content tables built by PNEUMA (see below), in case you want them outside the report.

</details>

[MultiQC](http://multiqc.info) summarises the whole run, not only read QC. Sections follow the pipeline flow:

| Order | Section | Source | Contents |
| --- | --- | --- | --- |
| 1 | General Statistics | all stages | one row per sample; PNEUMA adds `Reads after QC` (reads left after nanoq and Hostile, i.e. the input to taxonomy), `% Host removed` and `% Classified` (Kraken2) |
| 2 | NanoStat | NanoPlot `<sample>.NanoStats.txt` | read counts, N50, quality, length and quality distributions of the host-free reads |
| 3 | SeqKit stats | `preprocess/seqkit/` | sequence counts and length/quality statistics of the host-free reads |
| 4 | Porechop | `porechop_abi` logs | adapter trimming; only with `--run_porechop` |
| 5 | nanoq | `preprocess/nanoq/` | read-length and quality summary and distributions after filtering |
| 6 | Hostile host removal | `*.hostile.json` (custom table) | reads in, reads out, host reads, `% Host removed` per sample (MultiQC has no Hostile module) |
| 7 | Kraken2 | Kraken2 reports | classified/unclassified reads and top taxa |
| 8 | Bracken | Bracken reports | top taxa after re-estimation (the plot title says Kraken2 because Bracken writes the same report format) |
| 9 | Recentrifuge | `pneuma_controls.rcf.stat.csv` (custom table) | per-sample read and score statistics of the control-mode run; only when the samplesheet marks negative controls |
| 10 | ABRicate hits | `*_abricate.csv` (custom table) | hits per sample for `card`, `vfdb` and `plasmidfinder`, total hits and distinct genes; samples without hits show zeros |
| 11 | Software versions, workflow summary, methods | pipeline | versions of every tool and the parameters of the run |

The custom tables are written by [`bin/make_mqc_tables.py`](../bin/make_mqc_tables.py) (process `MULTIQC_TABLES`) and the section order is set in [`assets/multiqc_config.yml`](../assets/multiqc_config.yml).

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow, each suffixed with the launch timestamp (`<timestamp>` is `yyyy-MM-dd_HH-mm-ss`, from `params.trace_report_suffix`): `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt` and `pipeline_dag_<timestamp>.html`.
  - `params_<timestamp>.json`: the parameters of the run, written by the nf-core completion step.
  - `pneuma_software_mqc_versions.yml`: software versions of every process.
  - `database_manifest.tsv`: name, version, path and sha256 of each reference database (only when `--database_manifest` is given or `--db_dir` contains one). With `--build_databases` the manifest and the databases are written to `--db_outdir` instead.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt`; present only if `--email` / `--email_on_fail` is used.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) generates reports useful for troubleshooting: launch commands, run times and resource usage.

## Planned outputs

Phase 2 (priority-pathogen reconstruction and the WHO-priority respiratory panel) will add a pathogen report directory; it does not exist yet.
