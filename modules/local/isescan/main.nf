//
// ISEScan (insertion sequences, transposase HMMs; no database to stage) on the contigs that carry a hit.
//
process ISESCAN {
    tag 'hit_contigs'
    label 'process_medium'
    container 'quay.io/biocontainers/isescan:1.7.3--h7b50bb2_0'

    input:
    path fasta

    output:
    path 'hit_contigs.isescan.tsv', emit: tsv
    tuple val("${task.process}"), val('isescan'), eval("isescan.py --version 2>&1 | sed 's/.*isescan //I; s/ .*//'"), topic: versions, emit: versions_isescan

    script:
    def header = ['seqID', 'family', 'cluster', 'isBegin', 'isEnd', 'isLen', 'ncopy4is', 'start1', 'end1', 'start2', 'end2',
                  'score', 'irId', 'irLen', 'nGaps', 'orfBegin', 'orfEnd', 'strand', 'orfLen', 'E-value', 'E-value4copy',
                  'type', 'ov', 'tir'].join('\t')
    // ISEScan 1.7.3 writes <output>/<seqfile path>.tsv: isescan/<staged fasta name>.tsv
    """
    if [ -s ${fasta} ]; then
        isescan.py --seqfile ${fasta} --output isescan --nthread ${task.cpus}
    fi
    if [ -s isescan/${fasta}.tsv ]; then cp isescan/${fasta}.tsv hit_contigs.isescan.tsv; else printf '${header}\n' > hit_contigs.isescan.tsv; fi
    """

    stub:
    """
    printf 'seqID\tfamily\tcluster\tisBegin\tisEnd\tisLen\tstrand\n' > hit_contigs.isescan.tsv
    """
}
