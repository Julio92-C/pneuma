process ABRICATE_READS {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'

    input:
    tuple val(meta), path(fasta)
    path abricate_db    // optional: directory with current databases; [] = bundled ones

    output:
    tuple val(meta), path("${meta.id}_abricate.csv"), emit: csv
    path "versions.yml"                              , emit: versions

    script:
    def dbs = params.abricate_dbs.tokenize(',')
    def datadir = abricate_db ? "--datadir ${abricate_db}" : ''
    """
    : > ${meta.id}_abricate.csv
    if [ -s ${fasta} ]; then
        first=1
        for db in ${dbs.join(' ')}; do
            abricate ${fasta} --db \$db ${datadir} \
                --minid ${params.abricate_minid} --mincov ${params.abricate_mincov} \
                --threads ${task.cpus} --csv > abricate_\$db.tmp.csv
            # each ABRicate call prints its own header: keep only the first one
            if [ \$first -eq 1 ]; then cat abricate_\$db.tmp.csv; else tail -n +2 abricate_\$db.tmp.csv; fi >> ${meta.id}_abricate.csv
            first=0
        done
        rm -f abricate_*.tmp.csv
    else
        # no tagged reads: still write the single ABRicate header line (header only = no hits)
        echo '#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE' >> ${meta.id}_abricate.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(abricate --version | sed 's/abricate //')
        abricate_dbs: "\$(abricate --list ${datadir} | awk 'NR>1 {printf "%s:%s ", \$1, \$4}')"
    END_VERSIONS
    """
}
