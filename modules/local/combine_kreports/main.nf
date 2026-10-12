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
    // a single staged file arrives as a bare Path (iterating a Path yields its segments), so wrap it
    def report_list = reports instanceof List ? reports : [reports]
    def name_list   = names instanceof List ? names : [names]
    """
    combine_kreports.py \
        -r ${report_list.join(' ')} \
        -o ${prefix}.txt \
        --sample-names ${name_list.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: 1.2
    END_VERSIONS
    """
}
