process COMBINE_KREPORTS {
    tag "$prefix"
    label 'process_single'
    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'

    input:
    val  prefix
    val  names      // List<String>, same order as reports
    path reports, stageAs: 'reports/*'

    output:
    path "${prefix}.txt", emit: txt
    path "versions.yml" , emit: versions

    script:
    """
    combine_kreports.py \
        -r ${reports.join(' ')} \
        -o ${prefix}.txt \
        --sample-names ${names.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: 1.2
    END_VERSIONS
    """
}
