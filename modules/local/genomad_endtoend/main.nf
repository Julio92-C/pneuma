//
// geNomad end-to-end with a STAGED database, on one pooled FASTA (hit reads or contigs). Keeps the four tables the
// genomic-context analysis reads. An empty FASTA gives header-only tables (geNomad itself fails on empty input).
//
process GENOMAD_ENDTOEND {
    tag "$meta.id"
    label 'process_high'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.13.0--pyhdfd78af_0' :
        'quay.io/biocontainers/genomad:1.13.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${meta.id}_aggregated_classification.tsv"), emit: scores
    tuple val(meta), path("${meta.id}_virus_summary.tsv")            , emit: virus
    tuple val(meta), path("${meta.id}_plasmid_summary.tsv")          , emit: plasmids
    tuple val(meta), path("${meta.id}_genes.tsv")                    , emit: genes
    tuple val("${task.process}"), val('genomad'), eval("genomad --version 2>&1 | sed 's/.*version //'"), topic: versions, emit: versions_genomad

    script:
    def p = meta.id
    """
    if [ -s ${fasta} ]; then
        ln -sf ${fasta} ${p}.fna
        genomad end-to-end --threads ${task.cpus} --splits ${params.genomad_splits} ${p}.fna out ${db}
        cp out/${p}_aggregated_classification/${p}_aggregated_classification.tsv .
        cp out/${p}_summary/${p}_virus_summary.tsv out/${p}_summary/${p}_plasmid_summary.tsv .
        cp out/${p}_annotate/${p}_genes.tsv .
    else
        printf 'seq_name\\tchromosome_score\\tplasmid_score\\tvirus_score\\n' > ${p}_aggregated_classification.tsv
        printf 'seq_name\\tlength\\ttopology\\n' > ${p}_virus_summary.tsv
        printf 'seq_name\\tlength\\tn_hallmarks\\tconjugation_genes\\n' > ${p}_plasmid_summary.tsv
        printf 'gene\\tstart\\tend\\tplasmid_hallmark\\tannotation_conjscan\\tannotation_description\\n' > ${p}_genes.tsv
    fi
    """

    stub:
    def p = meta.id
    """
    printf 'seq_name\\tchromosome_score\\tplasmid_score\\tvirus_score\\n' > ${p}_aggregated_classification.tsv
    printf 'seq_name\\tlength\\ttopology\\n' > ${p}_virus_summary.tsv
    printf 'seq_name\\tlength\\tn_hallmarks\\tconjugation_genes\\n' > ${p}_plasmid_summary.tsv
    printf 'gene\\tstart\\tend\\tplasmid_hallmark\\tannotation_conjscan\\tannotation_description\\n' > ${p}_genes.tsv
    """
}
