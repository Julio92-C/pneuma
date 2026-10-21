process RECENTRIFUGE {
    label 'process_medium'
    container 'quay.io/biocontainers/recentrifuge:2.1.2--pyhdfd78af_0'

    input:
    val  controls     // List<String> sample ids of negative controls
    path assignments, stageAs: 'kraken/*'
    path taxdump

    output:
    path "pneuma.rcf.data.csv"   , emit: csv
    path "pneuma.rcf.html"       , emit: html
    path "pneuma_controls.rcf*.html", emit: controls_html, optional: true
    path "pneuma_controls.rcf.stat.csv", emit: stat, optional: true
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    ctrl_args=""
    other_args=""
    for f in kraken/*; do
        s=\$(basename "\$f"); s=\${s%%.*}
        case " ${controls.join(' ')} " in
            *" \$s "*) ctrl_args="\$ctrl_args -k \$f" ;;
            *)         other_args="\$other_args -k \$f" ;;
        esac
    done

    rcf -n ${taxdump} \$ctrl_args \$other_args -o pneuma.rcf.html ${args} -e CSV

    # Normalize Recentrifuge CSV headers for MERIDIAN v1.0.0 compatibility
    # Ensure row 1 col 1 = "Samples" and row 2 col 1 = "Stats" (Recentrifuge 2.1.2 leaves them empty)
    sed -i '1s/^,/Samples,/; 2s/^,/Stats,/' pneuma.rcf.data.csv

    if [ ${controls.size()} -gt 0 ]; then
        rcf -n ${taxdump} \$ctrl_args \$other_args -c ${controls.size()} -o pneuma_controls.rcf.html ${args} -e CSV

        # Normalize Recentrifuge CSV headers for control file as well
        sed -i '1s/^,/Samples,/; 2s/^,/Stats,/' pneuma_controls.rcf.data.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        recentrifuge: \$(rcf --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
