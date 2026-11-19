//
// Plasmid drafts grouped across samples into plasmid units with prevalence and abundance (bin/plasmid_units.py).
//
process PLASMID_UNITS {
    tag 'plasmids'
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path candidates
    path fasta
    path isescan
    path(recruit, stageAs: 'recruit/*')
    path taxdump          // directory with nodes.dmp + names.dmp, or [] (read taxa stay taxids)

    output:
    path 'plasmid_units.tsv'       , emit: units
    path 'plasmid_unit_members.tsv', emit: members
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    def tx = taxdump ? "--taxdump ${taxdump}" : ''
    """
    plasmid_units.py --candidates ${candidates} --fasta ${fasta} --isescan ${isescan} --recruit ${recruit} ${tx} \
        --min-identity ${params.plasmid_unit_min_identity} --out-prefix plasmid_units
    mv plasmid_units_members.tsv plasmid_unit_members.tsv
    """

    stub:
    """
    printf 'unit\tn_samples\tsamples\n' > plasmid_units.tsv
    printf 'unit\tcontig\tsample\n' > plasmid_unit_members.tsv
    """
}
