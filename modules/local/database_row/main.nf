process DATABASE_ROW {
    tag "$name"
    label 'process_single'
    label 'process_low_memory'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(name), val(order), val(version), val(rel), path(target)

    output:
    path "${order}_${name}.row.tsv", emit: row
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    def source = name == 'bakta_db' ? '--version-from bakta' :
                 name == 'amrfinderplus_db' ? '--version-from amrfinderplus' : "--version '${version}'"
    """
    database_row.py --name ${name} ${source} --rel '${rel}' --target ${target} --out ${order}_${name}.row.tsv
    """

    stub:
    """
    database_row.py --name ${name} --version stub --rel '${rel}' --target ${target} --out ${order}_${name}.row.tsv
    """
}
