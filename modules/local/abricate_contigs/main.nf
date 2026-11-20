//
// ABRicate on the pooled local-assembly contigs, one run per --context_databases database (same --abricate_db as
// the reads, so gene names match), plus contig lengths and the FASTA of the contigs carrying a hit.
// No contigs -> header-only hits table and empty FASTAs.
//
process ABRICATE_CONTIGS {
    tag 'contigs'
    label 'process_medium'
    container 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'   // same image as modules/local/abricate_reads

    input:
    path(contigs, stageAs: 'contigs/*')
    path abricate_db    // optional: directory with current databases; [] = bundled ones

    output:
    path 'all_contigs.fa'         , emit: fasta
    path 'all_contigs.lengths.tsv', emit: lengths
    path 'contig_hits.tsv'        , emit: hits
    path 'hit_contigs.fa'         , emit: hit_fasta
    tuple val("${task.process}"), val('abricate'), eval("abricate --version | sed 's/abricate //'"), topic: versions, emit: versions_abricate

    script:
    def dbs = params.context_databases.tokenize(',').join(' ')
    def datadir = abricate_db ? "--datadir ${abricate_db}" : ''
    """
    export HOME=\$PWD
    mkdir -p contigs
    find contigs/ -name '*.fa' | sort | xargs -r cat > all_contigs.fa
    awk '/^>/ {if (n) print n "\t" l; n = substr(\$1, 2); l = 0; next} {l += length(\$0)} END {if (n) print n "\t" l}' \
        all_contigs.fa > all_contigs.lengths.tsv
    printf '#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%%COVERAGE\t%%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n' > contig_hits.tsv
    if [ -s all_contigs.fa ]; then
        for db in ${dbs}; do
            abricate ${datadir} --db \$db --minid ${params.abricate_minid} --mincov ${params.abricate_mincov} \
                --threads ${task.cpus} --nopath --quiet all_contigs.fa > \$db.tab
            grep -v '^#' \$db.tab >> contig_hits.tsv || true
        done
    fi
    awk -F'\t' 'NR > 1 {print \$2}' contig_hits.tsv | sort -u > hit_contigs.ids
    awk 'FILENAME == ARGV[1] {keep[\$1]; next} /^>/ {p = (substr(\$1, 2) in keep)} p' hit_contigs.ids all_contigs.fa > hit_contigs.fa
    """

    stub:
    """
    : > all_contigs.fa; : > all_contigs.lengths.tsv; : > hit_contigs.fa
    printf '#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%%COVERAGE\t%%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n' > contig_hits.tsv
    """
}
