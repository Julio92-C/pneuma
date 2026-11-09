//
// The reads of one sample that carry an ABRicate hit, renamed "<sample>__<read>" so samples can be pooled for
// geNomad; the Kraken2 taxid moves to a side table (geNomad uses "|" in provirus names).
//
process HIT_READS {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(fasta), path(csv)

    output:
    tuple val(meta), path("${meta.id}.hit_reads.fa"), emit: fasta
    tuple val(meta), path("${meta.id}.hits.tsv")    , emit: hits
    tuple val(meta), path("${meta.id}.lengths.tsv") , emit: lengths
    tuple val(meta), path("${meta.id}.taxa.tsv")    , emit: taxa
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    hit_reads.py --sample ${meta.id} --fasta ${fasta} --abricate ${csv} --out-prefix ${meta.id}
    """

    stub:
    """
    hit_reads.py --sample ${meta.id} --fasta ${fasta} --abricate ${csv} --out-prefix ${meta.id}
    """
}
