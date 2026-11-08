process COVERAGE_GATE {
    tag "$meta.id"
    label 'process_single'
    label 'process_low_memory'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(coverage)   // samtools coverage table of the unit's reads vs the panel reference

    output:
    tuple val(meta), path("${meta.id}.gate.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    """
    coverage_gate.py --unit '${meta.id}' --coverage ${coverage} \\
        --min-depth ${params.gate_min_depth} --min-breadth ${params.gate_min_breadth} --out ${meta.id}.gate.tsv
    """
}
