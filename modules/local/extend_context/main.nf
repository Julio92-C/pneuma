//
// Final genomic context per read hit of one database: the call of the local-assembly contig at the gene's
// projected position when it is resolved, else the read call (bin/extend_context.py).
//
process EXTEND_CONTEXT {
    tag "$db"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(db), path(read_loci), path(contig_loci)
    path(pafs, stageAs: 'paf/*')

    output:
    tuple val(db), path("hit_context_${db}.tsv"), emit: context
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    extend_context.py --db ${db} --read-loci ${read_loci} --contig-loci ${contig_loci} --paf ${pafs} --out hit_context_${db}.tsv
    """

    stub:
    """
    printf 'sample\tseq\tdatabase\tgene\tfinal_context\n' > hit_context_${db}.tsv
    """
}
