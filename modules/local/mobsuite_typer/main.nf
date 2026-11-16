//
// MOB-suite mob_typer --multi on the hit contigs (each contig typed on its own: replicon, relaxase, mating-pair
// formation, predicted mobility) with a STAGED database directory (see MOBSUITE_RECON).
//
process MOBSUITE_TYPER {
    tag 'hit_contigs'
    label 'process_medium'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.9--pyhdfd78af_0' :
        'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0' }"

    input:
    path fasta
    path db      // mob_init database directory

    output:
    path 'hit_contigs.mobtyper.txt', emit: results
    tuple val("${task.process}"), val('mobsuite'), eval("mob_typer --version 2>&1 | sed 's/.*mob_typer //; s/ .*//'"), topic: versions, emit: versions_mobsuite

    script:
    """
    export HOME=\$PWD
    if [ -s ${fasta} ]; then
        mob_typer --multi --infile ${fasta} --out_file hit_contigs.mobtyper.txt --num_threads ${task.cpus} \
            --database_directory ${db}
    else
        printf 'sample_id\trep_type(s)\trelaxase_type(s)\tmpf_type\tpredicted_mobility\n' > hit_contigs.mobtyper.txt
    fi
    """

    stub:
    """
    printf 'sample_id\trep_type(s)\trelaxase_type(s)\tmpf_type\tpredicted_mobility\n' > hit_contigs.mobtyper.txt
    """
}
