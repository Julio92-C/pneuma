process HOSTILE_INDEX {
    tag "${fasta.name}"
    label 'process_medium'
    container 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'

    input:
    path fasta

    output:
    path "*.mmi"       , emit: mmi
    path "versions.yml", emit: versions

    script:
    def base = fasta.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
    """
    # Built once per run and shared by every HOSTILE_CLEAN task. Given a FASTA,
    # hostile builds this index itself next to the (staged) FASTA, i.e. once per
    # task in its work dir (~7.7 GB each for GRCh38). Same preset as hostile 2.0.2
    # (hostile/aligner.py: minimap2 -x map-ont -d <mmi> <fasta>).
    minimap2 -x map-ont -t ${task.cpus} -d ${base}.mmi ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def base = fasta.name.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
    """
    touch ${base}.mmi
    echo '"${task.process}": {minimap2: stub}' > versions.yml
    """
}
