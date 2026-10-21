process PREPARE_KRAKEN2_DB {
    tag "pluspf_${release}"
    label 'process_single'
    label 'process_low_memory'
    label 'process_long'
    label 'error_retry'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    val release   // YYYYMMDD of the Kraken2 PlusPF index (https://benlangmead.github.io/aws-indexes/k2)

    output:
    path "kraken2_pluspf_${release}", emit: db
    path "01_kraken2.row.tsv"        , emit: row
    path "versions.yml"              , emit: versions

    script:
    def url = "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_${release}.tar.gz"
    """
    mkdir -p kraken2_pluspf_${release}
    wget -nv -O k2_pluspf_${release}.tar.gz '${url}'
    sha=\$(sha256sum k2_pluspf_${release}.tar.gz | cut -d' ' -f1)
    [ -n "\$sha" ] || { echo "ERROR: empty sha256 for k2_pluspf_${release}.tar.gz" >&2; exit 1; }
    tar -xzf k2_pluspf_${release}.tar.gz -C kraken2_pluspf_${release}
    rm k2_pluspf_${release}.tar.gz    # the tarball is not published; its sha256 is kept in the manifest
    for f in hash.k2d opts.k2d taxo.k2d; do
        [ -s kraken2_pluspf_${release}/\$f ] || { echo "ERROR: \$f missing after extracting the Kraken2 archive" >&2; exit 1; }
    done
    printf 'kraken2_pluspf\\t%s\\t%s\\t%s\\n' '${release}' 'kraken2_pluspf_${release}' "\$sha" > 01_kraken2.row.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | sed -n '1s/^GNU Wget //; 1s/ .*//p')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p kraken2_pluspf_${release}
    for f in hash.k2d opts.k2d taxo.k2d; do echo stub > kraken2_pluspf_${release}/\$f; done
    sha=\$(echo stub | sha256sum | cut -d' ' -f1)
    printf 'kraken2_pluspf\\t%s\\t%s\\t%s\\n' '${release}' 'kraken2_pluspf_${release}' "\$sha" > 01_kraken2.row.tsv
    echo '"${task.process}": {wget: stub}' > versions.yml
    """
}
