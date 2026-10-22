process DATABASE_MANIFEST {
    label 'process_single'
    label 'process_low_memory'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path rows, stageAs: 'rows/*'          // one *.row.tsv per prepared database
    val  root                             // absolute db_outdir, prepended to every row's relative path
    path existing, stageAs: 'existing/*'  // previous database_manifest.tsv in db_outdir, or []

    output:
    path "database_manifest.tsv", emit: manifest
    path "versions.yml"         , emit: versions

    script:
    def keep = existing ? "--existing ${existing}" : ''
    """
    build_database_manifest.py --root '${root}' ${keep} --out database_manifest.tsv rows/*.row.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def keep_stub = existing ? "--existing ${existing}" : ''
    """
    build_database_manifest.py --root '${root}' ${keep_stub} --out database_manifest.tsv rows/*.row.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
