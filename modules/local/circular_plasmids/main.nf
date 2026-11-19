//
// Circular local-assembly contigs carrying a hit = draft complete plasmids (bin/circular_candidates.py).
//
process CIRCULAR_PLASMIDS {
    tag 'contigs'
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path(infos, stageAs: 'info/*')
    path hits
    path fasta
    path scores
    path plasmids
    path mobtyper      // or []
    path(pafs, stageAs: 'paf/*')
    path taxa          // or []

    output:
    path 'circular_plasmid_candidates.tsv', emit: tsv
    path 'circular_plasmid_candidates.fa' , emit: fasta
    path 'per_sample/*.circles.fa'        , emit: per_sample, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def mob = mobtyper ? "--mob-typer ${mobtyper}" : ''
    def tx = taxa ? "--taxa ${taxa}" : ''
    """
    circular_candidates.py --assembly-info ${infos} --hits ${hits} --fasta ${fasta} --genomad-scores ${scores} \
        --genomad-plasmids ${plasmids} ${mob} --paf ${pafs} ${tx} --out-prefix circular_plasmid_candidates
    """

    stub:
    """
    printf 'contig\\tsample\\tlength\\n' > circular_plasmid_candidates.tsv
    : > circular_plasmid_candidates.fa
    mkdir -p per_sample
    for f in ${infos}; do s=\$(basename "\$f"); s=\${s%%.*}; printf '>%s__contig_1\\nACGT\\n' "\$s" > "per_sample/\$s.circles.fa"; done
    """
}
