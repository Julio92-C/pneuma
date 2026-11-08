process RESOLVE_PANEL {
    tag "${panel.name}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path panel        // pathogen_panel.csv
    path references   // panel_references.tsv
    path taxdump      // directory with nodes.dmp, names.dmp (same release as the Kraken2 DB)

    output:
    path "panel_resolved.tsv", emit: tsv
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    """
    resolve_panel.py --panel ${panel} --references ${references} --taxdump ${taxdump} --out panel_resolved.tsv
    """

    stub:
    """
    printf 'species\\tslug\\ttier\\twho_category\\taction\\ttaxid\\tscientific_name\\tmatched_name_class\\ttaxids\\treference_accession\\tconfusable_with\\tconfusable_accessions\\tani_min\\tamrfinder_organism\\n' > panel_resolved.tsv
    """
}
