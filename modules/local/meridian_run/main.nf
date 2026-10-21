process MERIDIAN_RUN {
    tag "$group"
    label 'process_high'
    container { params.meridian_container }
    containerOptions { workflow.containerEngine == 'docker' ? '--entrypoint ""' : '' }

    input:
    tuple val(group), path(config)
    path metadata,         stageAs: 'inputs/*'
    path abricate_summary, stageAs: 'inputs/*'
    path kraken_combined,  stageAs: 'inputs/*'
    path bracken_arranged, stageAs: 'inputs/*'
    path rcf_csv,          stageAs: 'inputs/*'

    output:
    path "Datasets"     , emit: datasets
    path "Figures"      , emit: figures
    path "Reports"      , emit: reports
    // MERIDIAN >= 1.1.0 stage safety net: stages that failed but did not abort the run
    tuple val(group), path("Reports/stage_errors.csv"), optional: true, emit: stage_errors
    path "versions.yml" , emit: versions

    script:
    if (task.ext.test_noop) {
        // test-only (tests/meridian_noop_ok.config): no real MERIDIAN, just prove the wiring ran once per group
        return """
        mkdir -p Datasets Figures Reports
        touch Reports/pipeline_report.html
        echo "${group} ${config}" > Reports/group.txt
        ${task.ext.test_stage_error == group ? "printf 'stage,error\nnetwork,forced test failure\n' > Reports/stage_errors.csv" : ''}
        echo "\"${task.process}\": {meridian: noop}" > versions.yml
        """
    }
    """
    sed "s|__PROJECT_ROOT__|\$PWD|" ${config} > run_config.yaml
    mkdir -p Datasets Figures Reports
    # Configure renv to avoid sandbox activation in unwritable HOME (Singularity on HPC clusters)
    export RENV_CONFIG_SANDBOX_ENABLED=FALSE
    export RENV_PATHS_ROOT="\$PWD/.renv_root"
    unset R_LIBS R_LIBS_USER R_LIBS_SITE R_ENVIRON_USER
    # MERIDIAN renders templates/report.qmd IN PLACE and writes Rplots.pdf to the cwd, but /opt/meridian is
    # read-only under Singularity (and root-owned under docker). Build a writable mirror in the task dir:
    # symlinks to every top-level entry of /opt/meridian (the ~1 GB renv library stays shared and
    # read-only), except templates/ (real copy, Quarto writes there) and run_pipeline.R (real copy:
    # it derives repo_root from normalizePath(script), which would resolve a symlink back to /opt/meridian).
    RW="\$PWD/meridian_rw"
    mkdir -p "\$RW"
    for f in /opt/meridian/* /opt/meridian/.[!.]*; do
        b="\$(basename "\$f")"
        case "\$b" in
            templates|run_pipeline.R) cp -a "\$f" "\$RW/\$b" ;;
            *) ln -s "\$f" "\$RW/\$b" ;;
        esac
    done
    # Quarto/deno caches must not land in an unwritable HOME
    export XDG_CACHE_HOME="\$PWD/.xdg/cache" XDG_DATA_HOME="\$PWD/.xdg/data" XDG_CONFIG_HOME="\$PWD/.xdg/config"
    mkdir -p "\$XDG_CACHE_HOME" "\$XDG_DATA_HOME" "\$XDG_CONFIG_HOME"
    # renv: activate from the mirror's .Rprofile (renv/ is a symlink to /opt/meridian/renv)
    export R_PROFILE_USER="\$RW/.Rprofile"
    ( cd "\$RW" && Rscript -e 'lp <- .libPaths(); cat("libPaths:", lp, sep="\n"); if (!requireNamespace("ALDEx2", quietly=TRUE) || !any(grepl("/opt/meridian/renv/library|meridian_rw/renv/library", lp))) stop("renv did not activate /opt/meridian/renv/library")' )
    ( cd "\$RW" && Rscript run_pipeline.R "\$OLDPWD/run_config.yaml" )
    rm -rf "\$RW"

    # Quarto failures do not make run_pipeline.R exit non-zero: check the report ourselves
    if grep -Eq '^[[:space:]]+report:[[:space:]]+true' run_config.yaml && [ ! -s Reports/pipeline_report.html ]; then
        echo "ERROR: MERIDIAN stages.report is true but Reports/pipeline_report.html was not produced (Quarto render failed; see Reports/pipeline.log)" >&2
        exit 1
    fi

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
