//
// Genomic context of the hits of one ABRicate database (route reads or contigs): plasmid / chromosomal_mge /
// chromosome / unresolved. Writes the evidence-backed table and, for comparison, the geNomad score-only table.
//
process ARG_CONTEXT {
    tag "${route}:${db}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(route), val(db), path(hits), path(lengths), path(fasta), path(scores), path(virus), path(genes), path(mob)

    output:
    tuple val(route), val(db), path("loci_${route}_${db}.tsv")           , emit: loci
    tuple val(route), val(db), path("loci_${route}_${db}_score_only.tsv"), emit: score_only
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def mob_arg = mob ? (route == 'reads' ? "--mob-typer ${mob}" : "--mob-contigs ${mob}") : ''
    def flank = route == 'reads' ? params.context_min_flank : 0
    def common = "--route ${route} --abricate ${hits} --abricate-db ${db} --lengths ${lengths} --fasta ${fasta} " +
                 "--genomad-scores ${scores} --genomad-virus ${virus} --genomad-genes ${genes} ${mob_arg} " +
                 "--min-flank ${flank} --mge-window ${params.context_mge_window} " +
                 "--plasmid-min ${params.context_plasmid_min} --chrom-min ${params.context_chrom_min}"
    """
    arg_context.py ${common} --genomad-plasmids ${genes} --out loci_${route}_${db}.tsv
    arg_context.py ${common} --out loci_${route}_${db}_score_only.tsv
    """

    stub:
    """
    printf 'sample\\troute\\tseq\\tseq_len\\tgene\\tdrug_class\\tstart\\tend\\tcategory\\n' > loci_${route}_${db}.tsv
    cp loci_${route}_${db}.tsv loci_${route}_${db}_score_only.tsv
    """
}
