import json

from conftest import DATA, load_bin

mqc = load_bin("make_mqc_tables")

KRAKEN = (
    "25.00\t25\t25\tU\t0\tunclassified\n"
    "75.00\t75\t0\tR\t1\troot\n"
    "75.00\t75\t75\tD\t2\t  Bacteria\n"
)


def _hostile(path, reads_in, reads_out):
    path.write_text(json.dumps([{"version": "2.0.2", "reads_in": reads_in,
                                 "reads_out": reads_out, "reads_removed": reads_in - reads_out}]))
    return path


def test_parse_hostile_percent_removed(tmp_path):
    f = _hostile(tmp_path / "S1.hostile.json", 200, 150)
    data = mqc.parse_hostile([f])
    assert data["S1"]["reads_removed"] == 50
    assert data["S1"]["pct_host_removed"] == 25.0


def test_parse_hostile_skips_stub_and_zero_reads(tmp_path):
    stub = tmp_path / "STUB.hostile.json"
    stub.write_text("{}")
    zero = _hostile(tmp_path / "ZERO.hostile.json", 0, 0)
    data = mqc.parse_hostile([stub, zero])
    assert list(data) == ["ZERO"]
    assert data["ZERO"]["pct_host_removed"] == 0.0


def test_parse_kraken_percent_classified(tmp_path):
    f = tmp_path / "S1.kraken2.report.txt"
    f.write_text(KRAKEN)
    data = mqc.parse_kraken([f])
    assert data["S1"] == {"unclassified": 25, "classified": 75, "pct_classified": 75.0}


def test_parse_kraken_all_unclassified(tmp_path):
    f = tmp_path / "S2.kraken2.report.txt"
    f.write_text("100.00\t10\t10\tU\t0\tunclassified\n")
    assert mqc.parse_kraken([f])["S2"]["pct_classified"] == 0.0


def test_parse_abricate_counts_per_database(tmp_path):
    hits = tmp_path / "S1_abricate.csv"
    hits.write_text((DATA / "abricate_S1.csv").read_text())
    empty = tmp_path / "S2_abricate.csv"
    empty.write_text((DATA / "abricate_S2_empty.csv").read_text())
    data = mqc.parse_abricate([hits, empty])
    assert data["S2"]["total_hits"] == 0
    assert data["S1"]["total_hits"] == 2
    assert data["S1"]["card"] + data["S1"]["vfdb"] + data["S1"]["plasmidfinder"] == 2
    assert list(data["S2"])[:3] == ["card", "vfdb", "plasmidfinder"]


def test_parse_rcf_stat_transposes(tmp_path):
    f = tmp_path / "pneuma_controls.rcf.stat.csv"
    f.write_text(",kraken/A.kraken2.classifiedreads,kraken/B.kraken2.classifiedreads\n"
                 "Seqs. read,87.0,175.0\nSeqs. class.,80.0,170.0\n")
    data = mqc.parse_rcf_stat(f)
    assert data["A"]["Seqs. read"] == 87.0
    assert data["B"]["Seqs. class."] == 170.0


def test_main_writes_all_custom_content_files(tmp_path):
    h = _hostile(tmp_path / "S1.hostile.json", 100, 90)
    k = tmp_path / "S1.kraken2.report.txt"
    k.write_text(KRAKEN)
    a = tmp_path / "S1_abricate.csv"
    a.write_text((DATA / "abricate_S1.csv").read_text())
    r = tmp_path / "pneuma_controls.rcf.stat.csv"
    r.write_text(",kraken/S1.kraken2.classifiedreads\nSeqs. read,100.0\n")
    out = tmp_path / "out"
    mqc.main(["--hostile", str(h), "--kraken", str(k), "--abricate", str(a),
              "--rcf-stat", str(r), "--outdir", str(out)])
    assert sorted(p.name for p in out.iterdir()) == [
        "abricate_mqc.tsv", "hostile_mqc.tsv", "pneuma_general_stats_mqc.tsv", "recentrifuge_mqc.tsv"]

    general = (out / "pneuma_general_stats_mqc.tsv").read_text().splitlines()
    assert "# plot_type: 'generalstats'" in general
    table = [l for l in general if not l.startswith("#")]
    assert table == ["Sample\treads_after_qc\tpct_host_removed\tpct_classified",
                     "S1\t90\t10.00\t75.00"]

    hostile = (out / "hostile_mqc.tsv").read_text()
    assert "# id: 'hostile'" in hostile
    assert "S1\t100\t90\t10\t10.00" in hostile
    assert "# id: 'abricate_hits'" in (out / "abricate_mqc.tsv").read_text()


def test_main_without_inputs_writes_nothing(tmp_path):
    mqc.main(["--outdir", str(tmp_path / "o")])
    assert list((tmp_path / "o").iterdir()) == []


# --- failure visibility -------------------------------------------------------

def test_bad_hostile_json_warns_and_fails(tmp_path, capsys):
    bad = tmp_path / "S1.hostile.json"
    bad.write_text("{not json")
    rc = mqc.main(["--hostile", str(bad), "--outdir", str(tmp_path / "o")])
    err = capsys.readouterr().err
    assert "WARNING: %s not parsed: invalid JSON" % bad in err
    assert rc == 1


def test_hostile_without_reads_in_warns_and_fails(tmp_path, capsys):
    f = tmp_path / "S1.hostile.json"
    f.write_text('[{"version": "2.0.2"}]')
    rc = mqc.main(["--hostile", str(f), "--outdir", str(tmp_path / "o")])
    assert "WARNING: %s not parsed: no reads_in" % f in capsys.readouterr().err
    assert rc == 1


def test_empty_kraken_report_warns_and_fails(tmp_path, capsys):
    f = tmp_path / "S1.kraken2.report.txt"
    f.write_text("")
    rc = mqc.main(["--kraken", str(f), "--outdir", str(tmp_path / "o")])
    assert "WARNING: %s not parsed: empty" % f in capsys.readouterr().err
    assert rc == 1


def test_stub_hostile_json_is_not_an_error(tmp_path, capsys):
    stub = tmp_path / "S1.hostile.json"
    stub.write_text("{}")
    rc = mqc.main(["--hostile", str(stub), "--outdir", str(tmp_path / "o")])
    assert rc == 0
    assert "WARNING" not in capsys.readouterr().err


def test_partial_failure_warns_but_succeeds(tmp_path, capsys):
    good = _hostile(tmp_path / "S1.hostile.json", 100, 90)
    bad = tmp_path / "S2.hostile.json"
    bad.write_text("oops")
    rc = mqc.main(["--hostile", str(good), str(bad), "--outdir", str(tmp_path / "o")])
    assert rc == 0
    assert "S2.hostile.json not parsed" in capsys.readouterr().err
    assert (tmp_path / "o" / "hostile_mqc.tsv").exists()


# --- Recentrifuge sample names (real Myriad shape, anonymised IDs) -------------

RCF_STAT = (
    ",kraken/S1.kraken2.classifiedreads,kraken/S2.kraken2.classifiedreads,"
    "kraken/S3.kraken2.classifiedreads,kraken/CTRL1.kraken2.classifiedreads\n"
    "Seqs. read,21076.0,6204.0,4855.0,186528.0\n"
    "Seqs. unclass.,20000.0,6000.0,4800.0,186000.0\n"
    "Seqs. class.,1076.0,204.0,55.0,528.0\n"
    "Seqs. filtered,7344.0\n"
)


def test_parse_rcf_stat_real_shape_names_and_short_row(tmp_path):
    f = tmp_path / "pneuma_controls.rcf.stat.csv"
    f.write_text(RCF_STAT)
    data = mqc.parse_rcf_stat(f)
    assert list(data) == ["S1", "S2", "S3", "CTRL1"]
    assert data["CTRL1"]["Seqs. read"] == 186528.0
    assert data["S1"]["Seqs. filtered"] == 7344.0
    assert "Seqs. filtered" not in data["S2"]


def test_strip_directory_and_classifiedreads_suffixes():
    for raw in ("kraken/S1.kraken2.classifiedreads", "kraken/S1.kraken2.classifiedreads.txt",
                "S1.classifiedreads", "a/b/S1.classifiedreads.txt", "kraken\\S1.kraken2.classifiedreads"):
        assert mqc.sample_name(raw, mqc.RCF_SUFFIXES) == "S1"


def test_zero_byte_abricate_csv_warns_and_fails(tmp_path, capsys):
    f = tmp_path / "S1_abricate.csv"
    f.write_text("")
    rc = mqc.main(["--abricate", str(f), "--outdir", str(tmp_path / "o")])
    assert "WARNING: %s not parsed: empty file" % f in capsys.readouterr().err
    assert rc == 1
