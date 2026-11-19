//
// IS elements within --context_is_window bp of each ARG on a contig; composite transposons (bin/is_flanks.py).
//
process IS_FLANKS {
    tag 'contigs'
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path isescan
    path hits

    output:
    path 'is_flanks.tsv', emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    is_flanks.py --isescan ${isescan} --hits ${hits} --db card --window ${params.context_is_window} --out is_flanks.tsv
    """

    stub:
    """
    printf 'sample\tcontig\tgene\tcomposite_transposon\torientation\n' > is_flanks.tsv
    """
}
