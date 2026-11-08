//
// Medaka consensus with a STAGED model file. The nf-core medaka module has no model input, and the
// bioconda medaka 2.2.2 package does not bundle r1041_e82_400bps_bacterial_methylation, so it would
// download the model at run time: impossible on offline HPC compute nodes / read-only images.
//
process MEDAKA_CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:2.2.2--py312h3050eb1_0' :
        'quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0' }"

    input:
    tuple val(meta), path(reads), path(draft)
    path model        // <model>_model_pt.tar.gz (from --db_dir/medaka or --medaka_model_path)

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: assembly
    tuple val("${task.process}"), val('medaka'), eval('medaka --version 2>&1 | sed "s/medaka //"'), topic: versions, emit: versions_medaka

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export HOME=\$PWD    # torch / matplotlib caches; HOME is not writable for the container user
    medaka_consensus -t ${task.cpus} -m ${model} -i ${reads} -d ${draft} -o medaka_out ${args}
    mv medaka_out/consensus.fasta ${meta.id}.fasta
    """

    stub:
    """
    echo '>stub' > ${meta.id}.fasta
    echo 'ACGT' >> ${meta.id}.fasta
    """
}
