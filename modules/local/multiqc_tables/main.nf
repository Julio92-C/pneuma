process MULTIQC_TABLES {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path hostile_logs   // *.hostile.json
    path kraken_reports // *.kraken2.report.txt
    path abricate_csvs  // *_abricate.csv
    path rcf_stat       // pneuma_controls.rcf.stat.csv (optional: [] when there are no controls)
    path pathogen_cohort // pathogen_cohort.tsv (optional: [] without --pathogens)
    path stage_errors    // <comparison>.stage_errors.csv from MERIDIAN (optional: [])
    path hit_context     // context/hit_context.tsv (optional: [] with --genomic_context false)
    path plasmid_units   // mobilome/plasmid_units.tsv (optional: [] without the mobilome)
    path(context_loci, stageAs: 'context_loci/*')   // loci_reads_<db>{,_score_only}.tsv (optional: [])

    output:
    path "*_mqc.tsv"   , emit: tables
    path "versions.yml", emit: versions

    script:
    def hostile  = hostile_logs   ? "--hostile ${hostile_logs}"     : ''
    def kraken   = kraken_reports ? "--kraken ${kraken_reports}"    : ''
    def abricate = abricate_csvs  ? "--abricate ${abricate_csvs}"   : ''
    def rcf      = rcf_stat       ? "--rcf-stat ${rcf_stat}"        : ''
    def pathogen = pathogen_cohort ? "--pathogen-cohort ${pathogen_cohort}" : ''
    def errors   = stage_errors   ? "--meridian-stage-errors ${stage_errors}" : ''
    def context  = hit_context   ? "--hit-context ${hit_context}"     : ''
    def units    = plasmid_units ? "--plasmid-units ${plasmid_units}" : ''
    def loci     = context_loci  ? "--context-loci ${context_loci}"   : ''
    """
    make_mqc_tables.py ${hostile} ${kraken} ${abricate} ${rcf} ${pathogen} ${errors} ${context} ${units} ${loci} --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch hostile_mqc.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
