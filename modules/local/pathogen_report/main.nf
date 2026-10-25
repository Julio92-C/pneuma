process BUILD_PATHOGEN_REPORT {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path panel_resolved
    path candidates
    path samples_tsv
    path units_tsv
    path gate,          stageAs: 'gate/*'
    path species,       stageAs: 'species/*'
    path assemblies,    stageAs: 'assemblies/*'
    path mlst,          stageAs: 'mlst/*'
    path amr,           stageAs: 'amrfinderplus/*'
    path abricate,      stageAs: 'abricate/*'
    path mob,           stageAs: 'mobsuite/*'
    path mixed,         stageAs: 'mixed/*'
    path checkm2,       stageAs: 'checkm2/*'
    path quast,         stageAs: 'quast/*'
    path read_abricate, stageAs: 'functional/*'
    path panel_file
    path manifest       // database_manifest.tsv, or []

    output:
    path "pathogen_cohort.tsv", emit: cohort
    path "html/*.html"        , emit: html
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def manifest_arg = manifest ? "--manifest ${manifest}" : ''
    """
    pathogen_report.py \\
        --panel ${panel_resolved} --candidates ${candidates} --samples ${samples_tsv} --units ${units_tsv} \\
        --results . --panel-file ${panel_file} ${manifest_arg} \\
        --pipeline-version '${workflow.manifest.version}' \\
        --gate-min-depth ${params.gate_min_depth} --gate-min-breadth ${params.gate_min_breadth} --outdir .
    """
}
