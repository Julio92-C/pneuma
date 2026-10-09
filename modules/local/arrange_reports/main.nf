process ARRANGE_REPORTS {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path kraken_combined
    path bracken_combined

    output:
    path "kraken2_combined_reports.txt", emit: kraken
    path "bracken_arranged.csv"        , emit: bracken
    path "versions.yml"                , emit: versions

    script:
    """
    arrange_reports.py \
        --kraken ${kraken_combined} \
        --bracken ${bracken_combined} \
        --kraken-out kraken2_combined_reports.txt \
        --bracken-out bracken_arranged.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
