process MERIDIAN_CONFIG {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path template
    path metadata
    val  controls   // List<String>

    output:
    path "config.yaml" , emit: config
    path "versions.yml", emit: versions

    script:
    def re = params.meridian_random_effect ? "--random-effect ${params.meridian_random_effect}" : ''
    """
    make_meridian_config.py \
        --template ${template} --out config.yaml \
        --study-id ${params.meridian_study_id} \
        --metadata ${metadata} \
        --group ${params.meridian_group} ${re} \
        --controls "${controls.join(',')}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
