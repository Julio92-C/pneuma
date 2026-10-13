# PNEUMA usage

## 1. Stage databases (once)

Fetch a pinned release and run the script from the Nextflow assets directory (pin `-r <tag>` for reproducibility):

    nextflow pull Julio92-C/pneuma -r v0.1.1
    bash ~/.nextflow/assets/Julio92-C/pneuma/bin/prepare_databases.sh --outdir ~/Scratch/pneuma_db --kraken-release <YYYYMMDD> --taxdump-month <YYYY-MM>

The script downloads the pinned Kraken2 PlusPF index and NCBI taxdump, fetches the Hostile
`human-t2t-hla` index, and downloads **current** ABRicate databases (`card`, `vfdb`,
`plasmidfinder`) with `abricate-get_db` inside the same ABRicate container image the pipeline
uses (`quay.io/biocontainers/abricate:1.0.1--ha8f3691_1`), so database formats match. The
databases bundled in that image are dated 2021-Mar-27 and are used only when `--abricate_db`
is not set. It needs singularity (or apptainer, or docker) on the node.

It writes `database_manifest.tsv` (`name`, `version`, `path`, `sha256`) with one row per
database; for ABRicate the version is the database date and the sha256 is that of its
`sequences` file. Copy the five `--param` lines it prints (`--kraken2_db`, `--taxdump`,
`--host_index`, `--abricate_db`, `--database_manifest`). Use `--dry-run` to print the commands only. From a development clone use `bin/prepare_databases.sh` directly.

Use `~/Scratch/...` on Myriad: the home quota is too small for the Kraken2 index. Run
`prepare_databases.sh --help` for all options.

### Sizes and memory

| Database | Download | Disk after setup | Needed at run time |
| --- | --- | --- | --- |
| Kraken2 PlusPF | 85 GB (`.tar.gz`) | 110 GB (+ 85 GB tarball, deletable once extracted) | about 110 GB RAM |
| NCBI taxdump | 70 MB | under 1 GB | negligible |
| Hostile `human-t2t-hla` | 1 GB | a few GB | a few GB |
| ABRicate (3 DBs) | tens of MB | under 1 GB | negligible |

Sizes are approximate (Kraken2 figures from the Kraken2 index list for PlusPF, others from the
download URLs). Kraken2 loads the whole index into RAM, so it needs about the index size; the
`myriad` profile requests 110 GB for `KRAKEN2_KRAKEN2`, which matches PlusPF. Budget about 200 GB
of free disk while the script runs.

### Myriad job script for the download

The script needs little CPU or RAM but a lot of time and disk. Example (adjust the time, `-wd`
and the dates; `tmpfs` is the node-local scratch):

    #!/bin/bash -l
    #$ -N pneuma_db
    #$ -l h_rt=12:00:00
    #$ -l mem=4G
    #$ -l tmpfs=20G
    #$ -pe smp 1
    #$ -wd /home/<user>/Scratch/pneuma_db_logs
    module load singularity-env   # or whatever provides singularity on your system
    export NXF_SINGULARITY_CACHEDIR=~/Scratch/singularity_cache
    bash ~/.nextflow/assets/Julio92-C/pneuma/bin/prepare_databases.sh       --outdir ~/Scratch/pneuma_db --kraken-release <YYYYMMDD> --taxdump-month <YYYY-MM>

If compute nodes have no outbound internet, run the script on a login node inside `tmux` or
`screen` instead (so it survives a dropped connection), then use the printed parameters.

### Reusing an existing database

`--kraken2_db` takes a directory containing `hash.k2d`, `opts.k2d` and `taxo.k2d`. Bracken also
needs `database<READLEN>mers.kmer_distrib` in that directory, where `<READLEN>` matches
`--bracken_read_length` (default `300`); if your database only has other read lengths, set
`--bracken_read_length` to one of them (or build the file with `bracken-build`).

To fetch only the databases you do not have, skip the others. Skipped databases are not
downloaded, are left out of `database_manifest.tsv`, and no parameter line is printed for them;
pass your own path for them on the `nextflow run` command line. For example, with an existing
Kraken2 database:

    bash ~/.nextflow/assets/Julio92-C/pneuma/bin/prepare_databases.sh       --outdir ~/Scratch/pneuma_db --skip-kraken --taxdump-month <YYYY-MM>
    # then: --kraken2_db /path/to/existing/kraken2_db  plus the printed --taxdump, --host_index,
    #       --abricate_db and --database_manifest

(The manifest then lacks a Kraken2 row; add one by hand with the `name`, `version`, `path`,
`sha256` columns if you want it recorded.) The other switches are `--skip-taxdump`,
`--skip-hostile` and `--skip-abricate`; `--kraken-release` / `--taxdump-month` are only required
for databases that are not skipped.

## 2. Samplesheet

`sample,run_id,batch,site,fastq,is_control` — one row per sample × run. A sample sequenced on
several flow cells gets several rows with the same `sample`; their reads are concatenated.
`sample` must match `^[A-Za-z0-9][A-Za-z0-9-]*$`.

## 3. Run

    nextflow run Julio92-C/pneuma -r v0.1.1 -profile <docker|singularity|myriad> \
      --input samplesheet.csv --outdir results \
      --kraken2_db ... --taxdump ... --host_index ... --abricate_db ... --database_manifest ... \
      --meridian_metadata metadata.csv --meridian_group sample_type --meridian_random_effect ID \
      -resume

On UCL Myriad use `-profile myriad` (and `-profile test,myriad` for the test): it enables
Singularity and sets `singularity.cacheDir` to `~/Scratch/singularity_cache`. If you instead
combine `singularity` with another profile (for example `-profile test,singularity`), Nextflow warns
"Singularity cache directory has not been defined"; export `NXF_SINGULARITY_CACHEDIR=~/Scratch/singularity_cache`
before running to avoid repeated image pulls.

`--input` and `--outdir` are required; always pin `-r <tag>` for reproducibility. From a development
clone use `nextflow run . -profile ...`.

### ABRicate databases

`--abricate_db <dir>` (default `null`) points at the directory produced by
`prepare_databases.sh` (`<outdir>/abricate`, containing `card/`, `vfdb/`, `plasmidfinder/`).
When set, the directory is staged into each `ABRICATE_READS` task and passed as
`abricate --datadir`. When unset, the databases bundled in the container are used.
`--abricate_dbs` (default `card,vfdb,plasmidfinder`) is the comma-separated list of database *names* searched inside that directory.

## 4. Outputs

- `multiqc/` `multiqc_report.html` (QC and software versions)
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
