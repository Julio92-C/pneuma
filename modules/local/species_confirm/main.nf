process SPECIES_CONFIRM {
    tag "$meta.id"
    label 'process_single'
    label 'process_low_memory'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(skani)   // skani dist of the assembly vs panel reference + confusables
    path references                // panel_references.tsv (accession -> species)

    output:
    tuple val(meta), path("${meta.id}.species.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    """
    species_confirm.py --unit '${meta.id}' --species '${meta.species}' --ani-min ${meta.ani_min} \
        --skani ${skani} --references ${references} --out ${meta.id}.species.tsv
    """
}
