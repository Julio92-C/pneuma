process MERIDIAN_RUN {
    label 'process_high'
    container { params.meridian_container }
    containerOptions { workflow.containerEngine == 'docker' ? '--entrypoint ""' : '' }

    input:
    path config
    path metadata,         stageAs: 'inputs/*'
    path abricate_summary, stageAs: 'inputs/*'
    path kraken_combined,  stageAs: 'inputs/*'
    path bracken_arranged, stageAs: 'inputs/*'
    path rcf_csv,          stageAs: 'inputs/*'

    output:
    path "Datasets"     , emit: datasets
    path "Figures"      , emit: figures
    path "Reports"      , emit: reports
    path "versions.yml" , emit: versions

    script:
    """
    sed "s|__PROJECT_ROOT__|\$PWD|" ${config} > run_config.yaml
    mkdir -p Datasets Figures Reports
    # Configure renv to avoid sandbox activation in unwritable HOME (Singularity on HPC clusters)
    export RENV_CONFIG_SANDBOX_ENABLED=FALSE
    export RENV_PATHS_ROOT="\$PWD/.renv_root"
    export R_PROFILE_USER=/opt/meridian/.Rprofile
    unset R_LIBS R_LIBS_USER R_LIBS_SITE R_ENVIRON_USER
    ( cd /opt/meridian && Rscript run_pipeline.R "\$OLDPWD/run_config.yaml" )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meridian: ${params.meridian_container.tokenize(":").last()}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datasets Figures Reports
    touch Reports/pipeline_report.html
    echo "\"${task.process}\": {meridian: ${params.meridian_container.tokenize(":").last()}}" > versions.yml
    """
}
