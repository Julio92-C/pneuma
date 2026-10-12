process HOSTILE_CLEAN {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.hostile.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.hostile.json")    , emit: log
    path "versions.yml"                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # hostile creates a data dir under \$HOME even for custom indexes; HOME is read-only for the container user
    export HOME=\$PWD
    hostile clean \\
        --fastq1 ${reads} \\
        --aligner minimap2 \\
        --index ${index} \\
        --threads ${task.cpus} \\
        --output hostile_out \\
        --airplane \\
        ${args} > ${meta.id}.hostile.json
    mv hostile_out/*.fastq.gz ${meta.id}.hostile.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version 2>&1 | sed 's/^.*hostile //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo -n | gzip > ${meta.id}.hostile.fastq.gz
    echo '{}' > ${meta.id}.hostile.json
    echo '"${task.process}": {hostile: stub}' > versions.yml
    """
}
