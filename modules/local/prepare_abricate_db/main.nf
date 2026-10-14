process PREPARE_ABRICATE_DB {
    tag "card,vfdb,plasmidfinder"
    label 'process_single'
    label 'process_low_memory'
    label 'process_long'
    label 'error_retry'
    container 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'   // same image as modules/local/abricate_reads

    output:
    path "abricate"          , emit: db
    path "04_abricate.row.tsv", emit: row
    path "versions.yml"      , emit: versions

    script:
    """
    export HOME=\$PWD
    mkdir -p abricate
    # Current ABRicate databases, fetched with the image the pipeline runs (format match).
    for db in card vfdb plasmidfinder; do
        abricate-get_db --db \$db --dbdir abricate --force
    done
    # abricate-get_db only downloads the FASTA; --setupdb builds the BLAST indexes that
    # `abricate --datadir` needs at search time and registers the DB date.
    abricate --setupdb --datadir abricate
    abricate --list --datadir abricate > abricate_list.txt   # DATABASE SEQUENCES DBTYPE DATE

    : > 04_abricate.row.tsv
    for db in card vfdb plasmidfinder; do
        adate=\$(awk -v d=\$db '\$1==d {print \$4}' abricate_list.txt)
        [ -n "\$adate" ] || { echo "ERROR: could not parse the date of ABRicate database '\$db' from 'abricate --list'" >&2; exit 1; }
        asha=\$(sha256sum abricate/\$db/sequences | cut -d' ' -f1)
        [ -n "\$asha" ] || { echo "ERROR: empty sha256 for abricate/\$db/sequences" >&2; exit 1; }
        printf 'abricate_%s\\t%s\\tabricate/%s\\t%s\\n' "\$db" "\$adate" "\$db" "\$asha" >> 04_abricate.row.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(abricate --version | sed 's/abricate //')
    END_VERSIONS
    """

    stub:
    """
    for db in card vfdb plasmidfinder; do
        mkdir -p abricate/\$db
        echo ">stub_\$db" > abricate/\$db/sequences
    done
    : > 04_abricate.row.tsv
    for db in card vfdb plasmidfinder; do
        asha=\$(sha256sum abricate/\$db/sequences | cut -d' ' -f1)
        printf 'abricate_%s\\t2026-Jan-01\\tabricate/%s\\t%s\\n' "\$db" "\$db" "\$asha" >> 04_abricate.row.tsv
    done
    echo '"${task.process}": {abricate: stub}' > versions.yml
    """
}
