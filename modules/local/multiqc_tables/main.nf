process MULTIQC_TABLES {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path hostile_logs   // *.hostile.json
    path kraken_reports // *.kraken2.report.txt
    path abricate_csvs  // *_abricate.csv
    path rcf_stat       // pneuma_controls.rcf.stat.csv (optional: [] when there are no controls)

    output:
    path "*_mqc.tsv"   , emit: tables
    path "versions.yml", emit: versions

    script:
    def hostile  = hostile_logs   ? "--hostile ${hostile_logs}"     : ''
    def kraken   = kraken_reports ? "--kraken ${kraken_reports}"    : ''
    def abricate = abricate_csvs  ? "--abricate ${abricate_csvs}"   : ''
    def rcf      = rcf_stat       ? "--rcf-stat ${rcf_stat}"        : ''
    """
    make_mqc_tables.py ${hostile} ${kraken} ${abricate} ${rcf} --outdir .

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
