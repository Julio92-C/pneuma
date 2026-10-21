process RCF_FOR_MERIDIAN {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path rcf_csv
    val  samples   // List<String>: all samplesheet sample ids (controls included)

    output:
    path "pneuma.rcf.meridian.csv", emit: csv
    path "versions.yml"           , emit: versions

    script:
    """
    rcf_for_meridian.py --rcf ${rcf_csv} --samples "${samples.join(',')}" --out pneuma.rcf.meridian.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
