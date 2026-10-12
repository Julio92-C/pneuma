process ABRICATE_READS {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_abricate.csv"), emit: csv
    path "versions.yml"                              , emit: versions

    script:
    def dbs = params.abricate_dbs.tokenize(',')
    """
    : > ${meta.id}_abricate.csv
    if [ -s ${fasta} ]; then
        for db in ${dbs.join(' ')}; do
            abricate ${fasta} --db \$db \
                --minid ${params.abricate_minid} --mincov ${params.abricate_mincov} \
                --threads ${task.cpus} --csv >> ${meta.id}_abricate.csv
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(abricate --version | sed 's/abricate //')
        abricate_dbs: "\$(abricate --list | awk 'NR>1 {printf "%s:%s ", \$1, \$4}')"
    END_VERSIONS
    """
}
