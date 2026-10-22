# PNEUMA usage

## 1. Build the databases (once)

The primary route is the pipeline's own database mode. `--build_databases` together with
`--db_outdir` runs only the `PREPARE_DATABASES` subworkflow: no samplesheet is read and no analysis
is run. Pin `-r <tag>` for reproducibility:

    nextflow run Julio92-C/pneuma -r v0.2.0 -profile <docker|singularity|myriad> \
      --build_databases --db_outdir ~/Scratch/pneuma_db \
      --kraken_release <YYYYMMDD> --taxdump_month <YYYY-MM>

It downloads the pinned Kraken2 PlusPF index (`--kraken_release`, see
<https://benlangmead.github.io/aws-indexes/k2>) and NCBI taxdump (`--taxdump_month`, match the
Kraken2 release month), fetches the Hostile `human-t2t-hla` index, and downloads **current**
ABRicate databases (`card`, `vfdb`, `plasmidfinder`) with `abricate-get_db` inside the same ABRicate
container image the pipeline uses (`quay.io/biocontainers/abricate:1.0.1--ha8f3691_1`), so database
formats match. The databases bundled in that image are dated 2021-Mar-27 and are used only when no
ABRicate directory is supplied. Every download runs in a pinned container (Hostile
`2.0.2--pyhdfd78af_0`, ABRicate `1.0.1--ha8f3691_1`, Python `3.9--1` for the Kraken2/taxdump
downloads and the manifest), so the profile must provide Docker, Singularity or Apptainer.

Results are copied (`publishDir` mode `copy`) into `--db_outdir`:

    <db_outdir>/
      kraken2_pluspf_<YYYYMMDD>/      hash.k2d, opts.k2d, taxo.k2d (+ Bracken kmer_distrib files if shipped)
      taxdump_<YYYY-MM>/              nodes.dmp, names.dmp
      hostile/                        human-t2t-hla.mmi (+ the .fa.gz it was built from)
      abricate/                       card/, vfdb/, plasmidfinder/ (with BLAST indexes)
      database_manifest.tsv

`database_manifest.tsv` has the columns `name`, `version`, `path`, `sha256`, one row per database;
for ABRicate the version is the database date and the sha256 is that of its `sequences` file, for
Kraken2 and the taxdump the sha256 is that of the downloaded archive, and for Hostile that of the
`.mmi`. The run fails if any sha256 or date is empty. The Kraken2 tarball is deleted after extraction
(its checksum stays in the manifest).

Skip a database you already have with `--skip_kraken_download`, `--skip_taxdump_download`,
`--skip_hostile_download` or `--skip_abricate_download`; skipped databases are not downloaded and have
no manifest row. `--kraken_release` / `--taxdump_month` are only required for databases that are not
skipped. The download processes request little memory and a long time limit (`process_single`,
`process_low_memory`, `process_long`, with retries); the Hostile step builds a minimap2 index and uses
`process_low` (12 GB).

### Using the databases

    nextflow run Julio92-C/pneuma -r v0.2.0 -profile <docker|singularity|myriad> \
      --input samplesheet.csv --outdir results --db_dir ~/Scratch/pneuma_db ...

`--db_dir` fills `--kraken2_db`, `--taxdump`, `--host_index`, `--abricate_db` and
`--database_manifest` from that folder, using `database_manifest.tsv` (versions) and the layout
above. Parameters you pass explicitly override it, for example `--kraken2_db /path/to/existing/db`
or `--host_index my_host.fasta`. The required files are checked at start-up (`hash.k2d`, `opts.k2d`,
`taxo.k2d`; `nodes.dmp`, `names.dmp`; the `.mmi`; a `sequences` file for every database in
`--abricate_dbs`) and a clear error lists what is missing. `--kraken2_db`, `--taxdump` and
`--host_index` are required only when they cannot be resolved from `--db_dir`.

### Pathogen databases and `--skip_pathogens`

Pathogen reconstruction is opt-in: `--skip_pathogens` defaults to `true`, and a run without it needs no
pathogen database and behaves as in v0.2.0. With `--skip_pathogens false` the analysis run needs the
Phase 2 databases: `--checkm2_db`, `--bakta_db`, `--amrfinderplus_db`, `--medaka_model_path` and
`--panel_references_dir` (plus `--mobsuite_db` with `--run_mobsuite`). `--build_databases` builds them (CheckM2, Bakta
`--bakta_db_type light|full`, AMRFinderPlus, the Medaka model, an MLST scheme record, the panel
reference genomes and, unless `--skip_mobsuite_download`, MOB-suite) and `--db_dir` fills them. Use
`--skip_pathogen_download` to leave them out of a database build. With `--skip_pathogens false` and a
missing database the run stops with "--checkm2_db is required for pathogen reconstruction ...". `--pathogen_panel` and
`--panel_references` default to the files bundled in `assets/`.

### Fallback: `bin/prepare_databases.sh`

The same download as a stand-alone script, for systems where you would rather not run Nextflow for
it. Fetch it from the release tag you use (the Nextflow assets path is not stable: newer Nextflow
keeps per-commit clones under `~/.nextflow/assets/.repos/<org>/<repo>/clones/<sha>/`):

    curl -fsSL -o prepare_databases.sh https://raw.githubusercontent.com/Julio92-C/pneuma/v0.2.0/bin/prepare_databases.sh
    bash prepare_databases.sh --outdir ~/Scratch/pneuma_db --kraken-release <YYYYMMDD> --taxdump-month <YYYY-MM>

It writes the same layout and `database_manifest.tsv` (so `--db_dir ~/Scratch/pneuma_db` works for
its output too), prints the five `--param` lines, and needs singularity (or apptainer, or docker) on
the node. Use `--dry-run` to print the commands only; `prepare_databases.sh --help` lists all options
(`--skip-kraken`, `--skip-taxdump`, `--skip-hostile`, `--skip-abricate`). From a development clone use
`bin/prepare_databases.sh` directly.

Use `~/Scratch/...` on Myriad: the home quota is too small for the Kraken2 index.

### Command-line flags on Nextflow 26.04

Nextflow 26.04 enables the strict (v2) syntax parser by default, which passes command-line values as
strings, so boolean and numeric flags (`--build_databases`, `--skip_meridian`,
`--bracken_read_length 100`) fail schema validation with "Value is [string] but should be [boolean]".
This affects every flag of the pipeline, not only the database mode. Until the pipeline declares
typed parameters, either put such parameters in a `-params-file` (YAML/JSON keeps their types) or set
`NXF_SYNTAX_PARSER=v1` for the run. Nextflow 25.10 is not affected.

### Sizes and memory

| Database | Download | Disk after setup | Needed at run time |
| --- | --- | --- | --- |
| Kraken2 PlusPF | 85 GB (`.tar.gz`) | 110 GB (+ 85 GB tarball, deletable once extracted) | about 110 GB RAM |
| NCBI taxdump | 70 MB | under 1 GB | negligible |
| Hostile `human-t2t-hla` | 1 GB | a few GB | a few GB |
| ABRicate (3 DBs) | tens of MB | under 1 GB | negligible |

Sizes are approximate (Kraken2 figures from the Kraken2 index list for PlusPF, others from the
download URLs). Kraken2 loads the whole index into RAM, so it needs about the index size; the
`myriad` profile requests 110 GB for `KRAKEN2_KRAKEN2`, which matches PlusPF. Budget about 250 GB
of free disk while the databases are built (Nextflow's work directory holds the extracted
Kraken2 index and the published copy goes to `--db_outdir`; remove the work directory afterwards).

### Myriad job script for the download

The build needs little CPU or RAM but a lot of time and disk. Example (adjust the time, `-wd` and the
dates). Nextflow runs the processes on the node with Singularity, so keep the work directory on Scratch:

    #!/bin/bash -l
    #$ -N pneuma_db
    #$ -l h_rt=24:00:00
    #$ -l mem=8G
    #$ -pe smp 2
    #$ -wd /home/<user>/Scratch/pneuma_db_logs
    module load java   # or whatever provides Java and Nextflow on your system
    export NXF_SINGULARITY_CACHEDIR=~/Scratch/singularity_cache
    nextflow run Julio92-C/pneuma -r v0.2.0 -profile singularity \
      -work-dir ~/Scratch/pneuma_db_work \
      --build_databases --db_outdir ~/Scratch/pneuma_db \
      --kraken_release <YYYYMMDD> --taxdump_month <YYYY-MM>

(With Nextflow 26.04 put the `--` parameters in a `-params-file`, see above.) With `-profile myriad`
the download processes are submitted to SGE as separate jobs instead; use `-profile singularity` as
above when compute nodes have no outbound internet and run Nextflow inside `tmux` or `screen` on a
login node (so it survives a dropped connection), or use the fallback script there.

### Reusing an existing database

`--kraken2_db` takes a directory containing `hash.k2d`, `opts.k2d` and `taxo.k2d`. Bracken also
needs `database<READLEN>mers.kmer_distrib` in that directory, where `<READLEN>` matches
`--bracken_read_length` (default `300`); if your database only has other read lengths, set
`--bracken_read_length` to one of them (or build the file with `bracken-build`).

To fetch only the databases you do not have, skip the others. Skipped databases are not
downloaded and are left out of `database_manifest.tsv`; pass your own path for them. For example,
with an existing Kraken2 database:

    nextflow run Julio92-C/pneuma -r v0.2.0 -profile <docker|singularity|myriad> \
      --build_databases --db_outdir ~/Scratch/pneuma_db --skip_kraken_download --taxdump_month <YYYY-MM>
    # then, in analysis runs: --db_dir ~/Scratch/pneuma_db --kraken2_db /path/to/existing/kraken2_db

(The manifest then lacks a Kraken2 row; add one by hand with the `name`, `version`, `path`,
`sha256` columns if you want it recorded.) The fallback script has the matching `--skip-kraken`,
`--skip-taxdump`, `--skip-hostile` and `--skip-abricate` switches.

## 2. Samplesheet

`sample,run_id,batch,site,fastq,is_control` — one row per sample × run. A sample sequenced on
several flow cells gets several rows with the same `sample`; their reads are concatenated.
`sample` must match `^[A-Za-z0-9][A-Za-z0-9-]*$`.

## 3. Run

    nextflow run Julio92-C/pneuma -r v0.2.0 -profile <docker|singularity|myriad> \
      --input samplesheet.csv --outdir results \
      --db_dir ~/Scratch/pneuma_db \
      --meridian_metadata metadata.csv --meridian_group sample_type --meridian_random_effect ID \
      -resume

On UCL Myriad use `-profile myriad` (and `-profile test,myriad` for the test): it enables
Singularity and sets `singularity.cacheDir` to `~/Scratch/singularity_cache`. If you instead
combine `singularity` with another profile (for example `-profile test,singularity`), Nextflow warns
"Singularity cache directory has not been defined"; export `NXF_SINGULARITY_CACHEDIR=~/Scratch/singularity_cache`
before running to avoid repeated image pulls.

To compare more than one metadata column, give `--meridian_group` a comma-separated list, for
example `--meridian_group sample_type,Macrolide_prophylaxis`. MERIDIAN then runs once per column, in
parallel, with study id `<--meridian_study_id>_<column>` and results in `results/meridian/<column>/`.
Every listed column must exist in the metadata and have at least 2 distinct non-NA levels among the
non-control samples; otherwise the run stops at the start with an error naming the column.

`--input` and `--outdir` are required for an analysis run; always pin `-r <tag>` for reproducibility. From a development
clone use `nextflow run . -profile ...`.

### ABRicate databases

`--abricate_db <dir>` (default `null`) points at the directory produced by
`--build_databases` (`<db_outdir>/abricate`) or `prepare_databases.sh` (`<outdir>/abricate`), containing `card/`, `vfdb/`, `plasmidfinder/`).
When set (or taken from `--db_dir`), the directory is staged into each `ABRICATE_READS` task and passed as
`abricate --datadir`. When unset, the databases bundled in the container are used.
`--abricate_dbs` (default `card,vfdb,plasmidfinder`) is the comma-separated list of database *names* searched inside that directory.

## 4. Outputs

- `multiqc/` `multiqc_report.html` (QC, host removal, taxonomy, contamination, ABRicate hits and software versions)
- `preprocess/` nanoq stats, host-free reads, NanoPlot, seqkit stats
- `taxonomy/` Kraken2 reports + per-read assignments, Bracken, `bracken_arranged.csv`, `kraken2_combined_reports.txt`, Krona
- `contamination/` `pneuma.rcf.data.csv`, `pneuma.rcf.html`, `pneuma_controls.rcf.html`
- `functional/` per-sample ABRicate, `summary_reportClean.csv`
- `meridian/` MERIDIAN Datasets/Figures/Reports
- `pipeline_info/` timestamped Nextflow reports, `pneuma_software_mqc_versions.yml`, `database_manifest.tsv`

## Known Myriad issues

None recorded yet. The SGE headers generated by the `myriad` profile have to be checked on a
Myriad login node (`grep '^#\$' work/*/*/.command.run | sort -u`); if `qsub` rejects a resource
Nextflow adds itself, record the exact header line here.
