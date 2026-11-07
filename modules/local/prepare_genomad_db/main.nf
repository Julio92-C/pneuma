process PREPARE_GENOMAD_DB {
    tag "genomad 1.13.0"
    label 'process_low'
    label 'process_long'
    label 'error_retry'
    container 'quay.io/biocontainers/genomad:1.13.0--pyhdfd78af_0'   // same image as GENOMAD_ENDTOEND

    output:
    path "genomad_db"            , emit: dir
    path "14_genomad_db.row.tsv" , emit: row
    tuple val("${task.process}"), val('genomad'), eval("genomad --version 2>&1 | sed 's/.*version //'"), topic: versions, emit: versions_genomad

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    genomad download-database .
    ver=\$(cat genomad_db/version.txt 2>/dev/null || echo unknown)
    database_row.py --name genomad_db --version "\$ver" --rel genomad_db --target genomad_db --out 14_genomad_db.row.tsv
    """

    stub:
    """
    mkdir -p genomad_db
    echo stub > genomad_db/version.txt
    database_row.py --name genomad_db --version stub --rel genomad_db --target genomad_db --out 14_genomad_db.row.tsv
    """
}
