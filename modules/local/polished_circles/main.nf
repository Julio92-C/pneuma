//
// The circular plasmid drafts cut out of a sample's Medaka-polished contigs. Polishing against all of the sample's
// contigs keeps reads of a chromosomal copy of an IS element or gene on the chromosome instead of on the plasmid.
//
process POLISHED_CIRCLES {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'   // same image as modules/local/select_reads

    input:
    tuple val(meta), path(polished), path(circles)

    output:
    tuple val(meta), path("${meta.id}.circles.polished.fasta"), emit: fasta
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    script:
    """
    seqkit seq -n -i ${circles} > circles.ids
    seqkit grep -f circles.ids ${polished} > ${meta.id}.circles.polished.fasta
    want=\$(wc -l < circles.ids)
    got=\$(grep -c '^>' ${meta.id}.circles.polished.fasta || true)
    if [ "\$got" -ne "\$want" ]; then
        echo "WARNING: ${meta.id}: \$got of \$want circular drafts found in the polished contigs" >&2
    fi
    """

    stub:
    """
    cp ${circles} ${meta.id}.circles.polished.fasta
    """
}
