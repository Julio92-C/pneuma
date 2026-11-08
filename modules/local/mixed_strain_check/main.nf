process MIXED_STRAIN_CHECK {
    tag "$meta.id"
    label 'process_single'
    label 'process_low_memory'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(mpileup)   // pooled reads vs the pooled (patient-level) assembly

    output:
    tuple val(meta), path("${meta.id}.mixed.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    """
    mixed_strain_check.py --unit '${meta.id}' --mpileup ${mpileup} \
        --min-depth ${params.mixed_min_depth} --min-maf ${params.mixed_min_maf} --min-alt ${params.mixed_min_alt} \
        --max-het-per-kb ${params.mixed_max_het_per_kb} --min-evaluable ${params.mixed_min_evaluable} \
        --out ${meta.id}.mixed.tsv
    """
}
