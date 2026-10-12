process TAG_READS {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(classified)

    output:
    tuple val(meta), path("${meta.id}_kraken2ID_PFP.fasta"), emit: fasta
    path "versions.yml"                                     , emit: versions

    script:
    """
    tag_reads_taxid.py ${classified} ${meta.id}_kraken2ID_PFP.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
