process DATABASE_MANIFEST {
    label 'process_single'
    label 'process_low_memory'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path rows, stageAs: 'rows/*'   // one *.row.tsv per prepared database
    val  root                      // absolute db_outdir, prepended to every row's relative path

    output:
    path "database_manifest.tsv", emit: manifest
    path "versions.yml"         , emit: versions

    script:
    """
    build_database_manifest.py --root '${root}' --out database_manifest.tsv rows/*.row.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    build_database_manifest.py --root '${root}' --out database_manifest.tsv rows/*.row.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
