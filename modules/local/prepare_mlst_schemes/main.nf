process PREPARE_MLST_SCHEMES {
    tag "mlst 2.25.0"
    label 'process_single'
    label 'process_low_memory'
    container 'quay.io/biocontainers/mlst:2.25.0--hdfd78af_0'   // same image as nf-core MLST

    output:
    path "mlst"                    , emit: dir
    path "12_mlst_schemes.row.tsv" , emit: row
    tuple val("${task.process}"), val('mlst'), eval('mlst --version 2>&1 | sed "s/mlst //"'), topic: versions, emit: versions_mlst

    script:
    // The PubMLST schemes are bundled in the container (not downloaded): record which ones.
    """
    mkdir -p mlst
    fa="\$(mlst --help 2>&1 | sed -n "s/.*--blastdb.*default '\\([^']*\\)'.*/\\1/p")"
    [ -s "\$fa.nsq" ] || { echo "ERROR: bundled MLST BLAST database not found at \$fa" >&2; exit 1; }
    ver=\$(mlst --version 2>&1 | sed 's/^mlst //')
    {
        echo "# MLST schemes bundled in quay.io/biocontainers/mlst:2.25.0--hdfd78af_0 (mlst \$ver)"
        echo "# BLAST database (mlst.fa.nsq) sha256: \$(sha256sum "\$fa.nsq" | cut -d' ' -f1)"
        mlst --longlist
    } > mlst/mlst_schemes.txt
    sha=\$(sha256sum mlst/mlst_schemes.txt | cut -d' ' -f1)
    printf 'mlst_schemes\\t%s\\t%s\\t%s\\n' "\$ver" 'mlst/mlst_schemes.txt' "\$sha" > 12_mlst_schemes.row.tsv
    """

    stub:
    """
    mkdir -p mlst
    echo stub > mlst/mlst_schemes.txt
    sha=\$(sha256sum mlst/mlst_schemes.txt | cut -d' ' -f1)
    printf 'mlst_schemes\\tstub\\tmlst/mlst_schemes.txt\\t%s\\n' "\$sha" > 12_mlst_schemes.row.tsv
    """
}
