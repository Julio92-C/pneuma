process ABRICATE_MERGE {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path csvs, stageAs: 'abricate/*'

    output:
    path "summary_reportClean.csv", emit: summary
    path "versions.yml"           , emit: versions

    script:
    """
    merge_abricate.py --out summary_reportClean.csv abricate/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
