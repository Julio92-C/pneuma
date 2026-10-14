process PREPARE_TAXDUMP {
    tag "taxdump_${month}"
    label 'process_single'
    label 'process_low_memory'
    label 'process_long'
    label 'error_retry'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    val month   // YYYY-MM of the NCBI taxdump archive (match the Kraken2 release month)

    output:
    path "taxdump_${month}"  , emit: taxdump
    path "02_taxdump.row.tsv", emit: row
    path "versions.yml"      , emit: versions

    script:
    def url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${month}-01.zip"
    """
    mkdir -p taxdump_${month}
    wget -nv -O taxdmp_${month}-01.zip '${url}'
    sha=\$(sha256sum taxdmp_${month}-01.zip | cut -d' ' -f1)
    [ -n "\$sha" ] || { echo "ERROR: empty sha256 for taxdmp_${month}-01.zip" >&2; exit 1; }
    unzip -o -q taxdmp_${month}-01.zip nodes.dmp names.dmp -d taxdump_${month}
    printf 'ncbi_taxdump\\t%s\\t%s\\t%s\\n' '${month}' 'taxdump_${month}' "\$sha" > 02_taxdump.row.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | sed -n '1s/^GNU Wget //; 1s/ .*//p')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p taxdump_${month}
    echo stub > taxdump_${month}/nodes.dmp
    echo stub > taxdump_${month}/names.dmp
    sha=\$(echo stub | sha256sum | cut -d' ' -f1)
    printf 'ncbi_taxdump\\t%s\\t%s\\t%s\\n' '${month}' 'taxdump_${month}' "\$sha" > 02_taxdump.row.tsv
    echo '"${task.process}": {wget: stub}' > versions.yml
    """
}
