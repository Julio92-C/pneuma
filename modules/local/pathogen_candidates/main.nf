process PATHOGEN_CANDIDATES {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path panel_resolved
    path samples_tsv
    path kreports, stageAs: 'kreports/*'
    path bracken,  stageAs: 'bracken/*'     // [] when no sample has a Bracken table
    path seqkit,   stageAs: 'seqkit/*'      // seqkit stats tables (mean read length); [] when absent
    path refdir                             // panel reference genomes (<accession>.fna.gz): genome sizes

    output:
    path "pathogen_candidates.tsv", emit: candidates
    path "pathogen_samples.tsv"   , emit: samples
    tuple val("${task.process}"), val('python'), eval('python --version 2>&1 | sed "s/Python //"'), topic: versions, emit: versions_python

    script:
    // a single staged file arrives as a Path, not a List: normalise before joining
    def kr = kreports instanceof List ? kreports : [kreports]
    def br = bracken instanceof List ? bracken : [bracken]
    def bracken_arg = bracken ? "--bracken ${br.join(' ')}" : ''
    def sk = seqkit instanceof List ? seqkit : [seqkit]
    def seqkit_arg = seqkit ? "--seqkit-stats ${sk.join(' ')}" : ''
    """
    pathogen_candidates.py \\
        --panel ${panel_resolved} \\
        --samples ${samples_tsv} \\
        --kreports ${kr.join(' ')} \\
        ${bracken_arg} \\
        --min-reads ${params.pathogen_min_reads} \\
        ${seqkit_arg} \\
        --refdir ${refdir} \\
        --min-est-depth ${params.pathogen_min_est_depth} \\
        --out pathogen_candidates.tsv \\
        --samples-out pathogen_samples.tsv
    """
}
