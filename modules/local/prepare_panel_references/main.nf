process PREPARE_PANEL_REFERENCES {
    tag "${references.name}"
    label 'process_single'
    label 'process_low_memory'
    label 'error_retry'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path references   // assets/panel_references.tsv (or --panel_references)

    output:
    path "panel_references"            , emit: dir
    path "13_panel_references.row.tsv" , emit: row
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    """
    fetch_panel_references.py --references ${references} --outdir panel_references
    database_row.py --name panel_references --version "tsv-\$(sha256sum ${references} | cut -c1-12)" \\
        --rel panel_references --target panel_references/panel_references_index.tsv --out 13_panel_references.row.tsv
    """

    stub:
    """
    mkdir -p panel_references
    printf 'species\\taccession\\trole\\tfile\\tsha256\\n' > panel_references/panel_references_index.tsv
    database_row.py --name panel_references --version stub --rel panel_references \\
        --target panel_references/panel_references_index.tsv --out 13_panel_references.row.tsv
    """
}
