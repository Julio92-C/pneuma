//
// context/hit_context.tsv: every read hit of every context database with its final call, the read's Kraken2 taxon and
// the contig host vote (the MERIDIAN hand-off), and the context summaries (bin/taxon_gene_context.py,
// bin/arg_context_summary.py). Group summaries use --meridian_metadata and --meridian_group when both are set.
//
process CONTEXT_SUMMARY {
    tag 'context'
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path(tables, stageAs: 'tables/*')   // hit_context_<db>.tsv (local assembly) or loci_reads_<db>.tsv (layer 1 only)
    path taxa                           // hit_reads.taxa.tsv, or []
    path taxdump                        // directory with nodes.dmp and names.dmp
    path metadata                       // --meridian_metadata, or []

    output:
    path 'hit_context.tsv', emit: hit_context
    path 'summary/*.tsv'  , emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def files   = tables instanceof List ? tables : [tables]
    def hc      = files.findAll { f -> f.toString().tokenize('/').last().startsWith('hit_context_') }
    def rl      = files.findAll { f -> f.toString().tokenize('/').last().startsWith('loci_reads_') }
    def groups  = params.meridian_group ? params.meridian_group.toString().split(',').collect { g -> g.trim() }.findAll { g -> g }.join(',') : ''
    def hcArg   = hc ? "--hit-context ${hc.join(' ')}" : ''
    def rlArg   = rl ? "--read-loci ${rl.join(' ')}" : ''
    def txArg   = taxa ? "--taxa ${taxa}" : ''
    def byGroup = metadata && groups ? "--metadata ${metadata} --group-cols '${groups}'" : ''
    """
    mkdir -p summary
    taxon_gene_context.py ${hcArg} ${rlArg} ${txArg} --taxdump ${taxdump} \
        --out hit_context.tsv --resolution summary/taxon_resolution.tsv
    arg_context_summary.py --hit-context hit_context.tsv ${byGroup} --outdir summary
    """

    stub:
    """
    mkdir -p summary
    printf 'sample\tseq\tdatabase\telement\tgene\tfinal_context\n' > hit_context.tsv
    for t in arg_context_by_sample arg_context_by_class gene_context taxon_resolution; do printf 'stub\n' > summary/\$t.tsv; done
    """
}
