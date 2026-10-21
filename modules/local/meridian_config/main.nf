process MERIDIAN_CONFIG {
    tag "$group"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path template
    path metadata
    val  controls   // List<String>
    val  group      // one metadata column per MERIDIAN comparison

    output:
    tuple val(group), path("config.yaml"), emit: config
    path "versions.yml", emit: versions

    script:
    def safe = group.replaceAll(/[^A-Za-z0-9_.-]/, '_')
    def re = params.meridian_random_effect ? "--random-effect ${params.meridian_random_effect}" : ''
    """
    make_meridian_config.py \
        --template ${template} --out config.yaml \
        --study-id ${params.meridian_study_id}_${safe} \
        --metadata ${metadata} \
        --group "${group}" ${re} \
        --controls "${controls.join(',')}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
