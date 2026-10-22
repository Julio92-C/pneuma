process PREPARE_MEDAKA_MODEL {
    tag "$model"
    label 'process_single'
    label 'process_low_memory'
    label 'error_retry'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    val model    // medaka model name, e.g. r1041_e82_400bps_bacterial_methylation

    output:
    path "medaka"                  , emit: dir
    path "11_medaka_model.row.tsv" , emit: row
    tuple val("${task.process}"), val('wget'), eval("wget --help 2>&1 | sed -n 's/^BusyBox v\\([0-9.]*\\).*/busybox-\\1/p; s/^GNU Wget \\([0-9.]*\\).*/\\1/p' | head -1"), topic: versions, emit: versions_wget

    script:
    // Models are versioned with medaka: keep v2.2.2 in sync with the MEDAKA_CONSENSUS container.
    def url = "https://media.githubusercontent.com/media/nanoporetech/medaka/v2.2.2/medaka/data/${model}_model_pt.tar.gz"
    """
    mkdir -p medaka
    wget -nv -O medaka/${model}_model_pt.tar.gz '${url}'
    tar -tzf medaka/${model}_model_pt.tar.gz > /dev/null || { echo "ERROR: ${url} is not a medaka model archive (wrong --medaka_model?)" >&2; exit 1; }
    database_row.py --name medaka_model --version '${model}@medaka-2.2.2' \\
        --rel medaka/${model}_model_pt.tar.gz --target medaka/${model}_model_pt.tar.gz --out 11_medaka_model.row.tsv
    """

    stub:
    """
    mkdir -p medaka
    echo stub > medaka/${model}_model_pt.tar.gz
    database_row.py --name medaka_model --version stub --rel medaka/${model}_model_pt.tar.gz \\
        --target medaka/${model}_model_pt.tar.gz --out 11_medaka_model.row.tsv
    """
}
