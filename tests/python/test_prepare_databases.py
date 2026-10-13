import os
import stat
import subprocess
from conftest import BIN


def run(*args, env=None):
    return subprocess.run(["bash", str(BIN / "prepare_databases.sh"), *args],
                          capture_output=True, text=True, env=env)


def test_dry_run_prints_pinned_urls(tmp_path):
    r = run("--outdir", str(tmp_path), "--kraken-release", "20250402",
            "--taxdump-month", "2025-04", "--dry-run")
    assert r.returncode == 0, r.stderr
    assert "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20250402.tar.gz" in r.stdout
    assert "taxdmp_2025-04-01.zip" in r.stdout
    assert "hostile fetch --aligner minimap2 --name human-t2t-hla" in r.stdout
    assert "hostile:2.0.2--pyhdfd78af_0" in r.stdout


def test_dry_run_fetches_current_abricate_dbs(tmp_path):
    r = run("--outdir", str(tmp_path), "--kraken-release", "20250402",
            "--taxdump-month", "2025-04", "--dry-run")
    assert r.returncode == 0, r.stderr
    assert "abricate:1.0.1--ha8f3691_1" in r.stdout
    for db in ("card", "vfdb", "plasmidfinder"):
        assert f"abricate-get_db --db {db} --dbdir {tmp_path}/abricate" in r.stdout
    assert f"abricate --setupdb --datadir {tmp_path}/abricate" in r.stdout


def test_missing_args_fail(tmp_path):
    r = run("--outdir", str(tmp_path))
    assert r.returncode != 0
    assert "--kraken-release" in r.stderr


def _write_exe(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def _mock_env(tmp_path):
    """wget/tar/unzip/singularity mocked; returns (env, fakebin)."""
    fakebin = tmp_path / "fakebin"
    fakebin.mkdir()
    _write_exe(fakebin / "wget", '#!/usr/bin/env bash\nwhile [[ $# -gt 0 ]]; do [[ "$1" == -O ]] && { echo data > "$2"; shift; }; shift; done\n')
    _write_exe(fakebin / "tar", "#!/usr/bin/env bash\nexit 0\n")
    _write_exe(fakebin / "unzip", "#!/usr/bin/env bash\nexit 0\n")
    # fake `singularity exec IMAGE CMD...`: emulate hostile and the abricate tools
    _write_exe(fakebin / "singularity", r'''#!/usr/bin/env bash
echo "RAW $@" >> "$MOCK_LOG"
shift 4   # exec --bind DIR IMAGE
echo "$@" >> "$MOCK_LOG"
case "$1" in
  hostile) echo idx > "$HOSTILE_CACHE_DIR/human-t2t-hla.mmi" ;;
  abricate-get_db)
    while [[ $# -gt 0 ]]; do case "$1" in --db) db="$2";; --dbdir) dir="$2";; esac; shift; done
    mkdir -p "$dir/$db"; echo ">seq_$db" > "$dir/$db/sequences" ;;
  abricate)
    if [[ "$2" == "--list" ]]; then
      printf 'DATABASE\tSEQUENCES\tDBTYPE\tDATE\n'
      if [[ -n "$MOCK_NO_DATE" ]]; then printf 'card\t100\tnucl\n'; else printf 'card\t100\tnucl\t2026-Oct-01\n'; fi
      printf 'vfdb\t50\tnucl\t2026-Oct-02\nplasmidfinder\t10\tnucl\t2026-Oct-03\n'
    fi ;;
esac
''')
    env = dict(os.environ, PATH=f"{fakebin}:{os.environ['PATH']}",
               PNEUMA_CONTAINER_RT="singularity", MOCK_LOG=str(tmp_path / "calls.log"))
    return env, fakebin


ARGS = ("--kraken-release", "20250402", "--taxdump-month", "2025-04")


def test_manifest_records_abricate_dbs_with_mocked_container(tmp_path):
    """Full (non-dry) run with wget/tar/unzip/singularity mocked."""
    env, _ = _mock_env(tmp_path)
    out = tmp_path / "db"
    r = run("--outdir", str(out), "--kraken-release", "20250402", "--taxdump-month", "2025-04", env=env)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in (out / "database_manifest.tsv").read_text().splitlines()]
    assert rows[0] == ["name", "version", "path", "sha256"]
    by_name = {row[0]: row for row in rows[1:]}
    assert set(by_name) == {"kraken2_pluspf", "ncbi_taxdump", "hostile_index",
                            "abricate_card", "abricate_vfdb", "abricate_plasmidfinder"}
    assert by_name["abricate_card"][1] == "2026-Oct-01"
    assert by_name["abricate_plasmidfinder"][1] == "2026-Oct-03"
    assert len(by_name["abricate_vfdb"][3]) == 64
    assert by_name["abricate_vfdb"][2] == str(out / "abricate" / "vfdb")
    assert f"--abricate_db {out}/abricate" in r.stdout
    assert "--database_manifest" in r.stdout
    calls = (tmp_path / "calls.log").read_text()
    assert "abricate --setupdb --datadir" in calls
    assert f"--bind {out}" in calls


def test_outdir_with_spaces_works(tmp_path):
    env, _ = _mock_env(tmp_path)
    out = tmp_path / "my db dir"
    r = run("--outdir", str(out), *ARGS, env=env)
    assert r.returncode == 0, r.stderr
    assert (out / "database_manifest.tsv").exists()
    calls = (tmp_path / "calls.log").read_text()
    assert "abricate --setupdb --datadir" in calls
    assert (out / "abricate" / "card" / "sequences").exists()


def test_empty_sha_fails(tmp_path):
    env, fakebin = _mock_env(tmp_path)
    _write_exe(fakebin / "sha256sum", "#!/usr/bin/env bash\nexit 0\n")   # prints nothing
    r = run("--outdir", str(tmp_path / "db"), *ARGS, env=env)
    assert r.returncode != 0
    assert "sha256" in r.stderr


def test_empty_abricate_date_fails(tmp_path):
    env, _ = _mock_env(tmp_path)
    env["MOCK_NO_DATE"] = "1"
    r = run("--outdir", str(tmp_path / "db"), *ARGS, env=env)
    assert r.returncode != 0
    assert "date" in r.stderr and "card" in r.stderr
    assert not (tmp_path / "db" / "database_manifest.tsv").exists() or         "unknown" not in (tmp_path / "db" / "database_manifest.tsv").read_text()


def test_missing_hostile_index_fails(tmp_path):
    env, fakebin = _mock_env(tmp_path)
    # hostile mock writes nothing -> find matches nothing
    mock = (fakebin / "singularity").read_text().replace(
        'hostile) echo idx > "$HOSTILE_CACHE_DIR/human-t2t-hla.mmi" ;;', 'hostile) : ;;')
    _write_exe(fakebin / "singularity", mock)
    r = run("--outdir", str(tmp_path / "db"), *ARGS, env=env)
    assert r.returncode != 0
    assert "human-t2t-hla" in r.stderr


def test_help_prints_usage_and_exits_zero():
    r = run("--help")
    assert r.returncode == 0
    assert "--skip-kraken" in r.stdout and "Usage" in r.stdout


def test_skip_kraken_and_taxdump_need_no_release_and_download_nothing(tmp_path):
    r = run("--outdir", str(tmp_path), "--skip-kraken", "--skip-taxdump", "--dry-run")
    assert r.returncode == 0, r.stderr
    assert "k2_pluspf" not in r.stdout and "taxdmp" not in r.stdout
    assert "hostile fetch" in r.stdout and "abricate-get_db" in r.stdout


def test_skip_all_dry_run_does_nothing_but_mkdir(tmp_path):
    r = run("--outdir", str(tmp_path), "--skip-kraken", "--skip-taxdump",
            "--skip-hostile", "--skip-abricate", "--dry-run")
    assert r.returncode == 0, r.stderr
    assert "wget" not in r.stdout and "hostile" not in r.stdout and "abricate" not in r.stdout


def test_skipped_databases_are_omitted_from_manifest_and_params(tmp_path):
    env, _ = _mock_env(tmp_path)
    out = tmp_path / "db"
    r = run("--outdir", str(out), "--skip-kraken", "--skip-taxdump", env=env)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in (out / "database_manifest.tsv").read_text().splitlines()]
    names = {row[0] for row in rows[1:]}
    assert names == {"hostile_index", "abricate_card", "abricate_vfdb", "abricate_plasmidfinder"}
    assert "--kraken2_db" not in r.stdout and "--taxdump" not in r.stdout
    assert "--host_index" in r.stdout and "--database_manifest" in r.stdout


def test_skip_hostile_and_abricate_keeps_kraken_and_taxdump(tmp_path):
    env, _ = _mock_env(tmp_path)
    out = tmp_path / "db"
    r = run("--outdir", str(out), *ARGS, "--skip-hostile", "--skip-abricate", env=env)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in (out / "database_manifest.tsv").read_text().splitlines()]
    assert {row[0] for row in rows[1:]} == {"kraken2_pluspf", "ncbi_taxdump"}
    assert "--host_index" not in r.stdout and "--abricate_db" not in r.stdout
    assert not (tmp_path / "calls.log").exists()
