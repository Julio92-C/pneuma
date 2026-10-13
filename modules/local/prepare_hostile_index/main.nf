process PREPARE_HOSTILE_INDEX {
    tag "human-t2t-hla"
    label 'process_low'
    label 'process_long'
    label 'error_retry'
    container 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'   // same tag as modules/local/hostile_clean

    output:
    path "hostile"           , emit: index
    path "03_hostile.row.tsv", emit: row
    path "versions.yml"      , emit: versions

    script:
    """
    # hostile writes under \$HOME even with a custom cache dir; HOME is read-only for the container user
    export HOME=\$PWD
    export HOSTILE_CACHE_DIR=\$PWD/hostile
    mkdir -p hostile
    hostile index fetch --name human-t2t-hla --minimap2

    mmi=\$(find hostile -name 'human-t2t-hla*.mmi' | sort | sed -n 1p)
    [ -n "\$mmi" ] || { echo "ERROR: no human-t2t-hla*.mmi found under hostile/ (hostile index fetch failed?)" >&2; exit 1; }
    sha=\$(sha256sum "\$mmi" | cut -d' ' -f1)
    [ -n "\$sha" ] || { echo "ERROR: empty sha256 for \$mmi" >&2; exit 1; }
    printf 'hostile_index\\thuman-t2t-hla\\t%s\\t%s\\n' "\$mmi" "\$sha" > 03_hostile.row.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version 2>&1 | sed 's/^.*hostile //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p hostile
    echo stub > hostile/human-t2t-hla.mmi
    sha=\$(echo stub | sha256sum | cut -d' ' -f1)
    printf 'hostile_index\\thuman-t2t-hla\\thostile/human-t2t-hla.mmi\\t%s\\n' "\$sha" > 03_hostile.row.tsv
    echo '"${task.process}": {hostile: stub}' > versions.yml
    """
}
