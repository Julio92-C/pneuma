//
// One sample's Flye contigs renamed "<sample>__contig_N", with an assembly-info table (length, coverage, circular).
//
process CONTEXT_CONTIGS {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(fasta), path(info)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa")       , emit: fasta
    tuple val(meta), path("${meta.id}.assembly_info.tsv"), emit: info
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    context_contigs.py --sample ${meta.id} --fasta ${fasta} --info ${info} --out-prefix ${meta.id}
    """

    stub:
    """
    : > ${meta.id}.contigs.fa
    printf 'sample\tcontig\tlength\tcoverage\tcircular\n' > ${meta.id}.assembly_info.tsv
    """
}
