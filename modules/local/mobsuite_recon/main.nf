//
// MOB-suite mob_recon with a STAGED database directory (the nf-core module cannot stage one; without
// it mob_recon downloads into its package directory, impossible on read-only/offline HPC nodes).
//
process MOBSUITE_RECON {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.9--pyhdfd78af_0' :
        'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db      // mob_init database directory

    output:
    tuple val(meta), path("${meta.id}.contig_report.txt")   , emit: contig_report
    tuple val(meta), path("${meta.id}.mobtyper_results.txt"), emit: mobtyper, optional: true
    tuple val("${task.process}"), val('mobsuite'), eval("mob_recon --version 2>&1 | sed 's/.*mob_recon //; s/ .*//'"), topic: versions, emit: versions_mobsuite

    script:
    def args = task.ext.args ?: ''
    """
    export HOME=\$PWD
    mob_recon --infile ${fasta} --outdir results --num_threads ${task.cpus} --sample_id ${meta.id} \
        --database_directory ${db} ${args}
    cp results/contig_report.txt ${meta.id}.contig_report.txt
    if [ -f results/mobtyper_results.txt ]; then cp results/mobtyper_results.txt ${meta.id}.mobtyper_results.txt; fi
    """
}
