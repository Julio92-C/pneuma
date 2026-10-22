process PREPARE_MOBSUITE_DB {
    tag "mob_suite 3.1.9"
    label 'process_low'
    label 'process_long'
    label 'error_retry'
    container 'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0'   // same image as MOBSUITE_RECON

    output:
    path "mobsuite"               , emit: dir
    path "10_mobsuite_db.row.tsv" , emit: row
    tuple val("${task.process}"), val('mobsuite'), eval("mob_init --version 2>&1 | sed 's/.*mob_init //; s/ .*//'"), topic: versions, emit: versions_mobsuite

    script:
    """
    export HOME=\$PWD
    mob_init -d mobsuite
    ver=\$(mob_init --version 2>&1 | sed 's/.*mob_init //; s/ .*//')
    database_row.py --name mobsuite_db --version "mob_suite-\$ver" --rel mobsuite --target mobsuite --out 10_mobsuite_db.row.tsv
    """

    stub:
    """
    mkdir -p mobsuite
    echo stub > mobsuite/stub.txt
    database_row.py --name mobsuite_db --version stub --rel mobsuite --target mobsuite --out 10_mobsuite_db.row.tsv
    """
}
