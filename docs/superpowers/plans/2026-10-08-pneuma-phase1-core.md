# PNEUMA Phase 1 (Core Pipeline) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A private Nextflow pipeline `Julio92-C/pneuma` that turns raw ONT FASTQ (multi-run, two cohorts) into QC'd host-free reads, Kraken2/Bracken profiles, Recentrifuge contamination output, taxid-tagged ABRicate calls and a MERIDIAN cohort analysis — replacing every script in `JC_2026-09_EBSP/Bash_scripts`.

**Architecture:** nf-core-template DSL2 pipeline. Five local subworkflows (PREPROCESS, TAXONOMY, CONTAMINATION, FUNCTIONAL_READS, MERIDIAN) wired in `workflows/pneuma.nf`. nf-core modules for standard tools; small local modules for glue; stdlib-only Python helpers in `bin/` (unit-tested with pytest); MERIDIAN runs from its published container.

**Tech Stack:** Nextflow ≥ 24.10.0, nf-core tools ≥ 3.3, nf-schema, nf-test ≥ 0.9, Python 3.9 (stdlib only), pytest, Docker (dev, WSL2) / Singularity (Myriad, SGE).

**Spec:** `docs/superpowers/specs/2026-10-08-pneuma-design.md` (copied into the pneuma repo in Task 1).

**Follow-on plans (written after this one ships):** Phase 2 — pathogen reconstruction + pathogen report (spec §4.2, §5.5, §5.6, §8.3). Phase 3 — regression vs. legacy outputs, spike-ins, full two-year cohort run (spec §8.2–8.4).

## Global Constraints

- Repo `Julio92-C/pneuma`, **private**, MIT license (`Copyright (c) 2026 Julio Cesar Ortega Cambara`).
- Nextflow ≥ 24.10.0; nf-core tools ≥ 3.3; nf-test ≥ 0.9.
- Development happens inside WSL2 (Ubuntu) — never run Nextflow from Windows PowerShell.
- Sample IDs match `^[A-Za-z0-9][A-Za-z0-9-]*$` (no `_`, no `.`).
- `bin/` helpers: Python 3.9-compatible, standard library only, executable (`chmod +x`), shebang `#!/usr/bin/env python3`.
- MERIDIAN container pinned to `ghcr.io/julio92-c/meridian:v1.0.0`.
- File-name contracts consumed by MERIDIAN v1.0.0 (do not rename): `<sample>_kraken2ID_PFP.fasta` (ABRicate `#FILE`), read headers `<read_id>_kraken:taxid|<taxid>`, `summary_reportClean.csv`, `bracken_arranged.csv`, `kraken2_combined_reports.txt`, `pneuma.rcf.data.csv`.
- ABRicate: databases `card,vfdb,plasmidfinder`, `--minid 80 --mincov 80`.
- Recentrifuge flags: `-s NORMA -y 25 -x 9606 -d -e CSV`.
- Kraken2 production DB: full PlusPF, one pinned release.
- Host removal: Hostile, minimap2 aligner, `human-t2t-hla` index.
- Every local module emits `versions.yml` in nf-core format.

---

## File map

```
pneuma/
├── LICENSE                                   MIT
├── main.nf / nextflow.config / nextflow_schema.json / modules.json   (template, edited)
├── workflows/pneuma.nf                       wires subworkflows
├── subworkflows/local/
│   ├── utils_nfcore_pneuma_pipeline/main.nf  (template) samplesheet parsing, grouping, control warning
│   ├── preprocess/main.nf                    cat → porechop_abi? → nanoq → hostile → nanoplot/seqkit, LOW_DEPTH
│   ├── taxonomy/main.nf                      kraken2 → bracken → combine_kreports → arrange → krona
│   ├── contamination/main.nf                 recentrifuge
│   ├── functional_reads/main.nf              tag_reads → abricate_reads → abricate_merge
│   └── meridian/main.nf                      meridian_config → meridian_run
├── modules/local/
│   ├── hostile_clean/main.nf
│   ├── combine_kreports/main.nf
│   ├── arrange_reports/main.nf
│   ├── recentrifuge/main.nf
│   ├── tag_reads/main.nf
│   ├── abricate_reads/main.nf
│   ├── abricate_merge/main.nf
│   ├── meridian_config/main.nf
│   └── meridian_run/main.nf
├── modules/nf-core/   cat/fastq, porechop/abi, nanoq, nanoplot, seqkit/stats, untar,
│                      kraken2/kraken2, bracken/bracken, krakentools/kreport2krona,
│                      krona/ktimporttext, multiqc
├── bin/
│   ├── tag_reads_taxid.py
│   ├── merge_abricate.py
│   ├── arrange_reports.py
│   ├── make_meridian_config.py
│   └── prepare_databases.sh
├── assets/
│   ├── schema_input.json
│   ├── samplesheet_test.csv
│   └── meridian_config.template.yaml
├── conf/ base.config · modules.config · test.config · myriad.config
├── tests/
│   ├── python/ conftest.py · test_tag_reads_taxid.py · test_merge_abricate.py
│   │           · test_arrange_reports.py · test_make_meridian_config.py · test_prepare_databases.py
│   ├── data/   classified_small.fastq · abricate_*.csv · bracken_combined.txt
│   │           · kraken_combined.txt · metadata_test.csv · taxdump/ (generated in Task 8)
│   ├── make_test_taxdump.py
│   └── main.nf.test
└── docs/ usage.md · output.md · superpowers/specs/2026-10-08-pneuma-design.md
```

---

### Task 1: Scaffold the repository

**Files:**
- Create: whole template tree, `LICENSE`, `docs/superpowers/specs/2026-10-08-pneuma-design.md`, `.github/workflows/pytest.yml`, `tests/python/conftest.py`

**Interfaces:**
- Produces: a template pipeline that runs `-profile test,docker`; pytest harness importing `bin/*.py` as modules.

- [ ] **Step 1: Install tooling in WSL2**

```bash
sudo apt-get update && sudo apt-get install -y openjdk-17-jre-headless python3-pip python3-venv git
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
python3 -m venv ~/.venvs/pneuma && source ~/.venvs/pneuma/bin/activate
pip install "nf-core>=3.3" pytest
curl -fsSL https://get.nf-test.com | bash && sudo mv nf-test /usr/local/bin/
nextflow -version && nf-core --version && nf-test version
```
Expected: Nextflow ≥ 24.10.0, nf-core ≥ 3.3, nf-test ≥ 0.9. Docker Desktop with WSL2 integration must be on (`docker run hello-world` works).

- [ ] **Step 2: Create the pipeline from the nf-core template**

Write `~/pneuma-template.yml`:
```yaml
name: pneuma
description: Long-read metagenomics of respiratory samples — QC, host removal, taxonomy, contamination, AMR/VF/MGE, priority-pathogen reconstruction and MERIDIAN cohort analysis
author: Julio Cesar Ortega Cambara
org: Julio92-C
outdir: .
is_nfcore: false
skip_features:
  - igenomes
  - fastqc
  - gitpod
  - codespaces
  - seqera_platform
  - adaptivecard
  - slackreport
  - rocrate
```
Run:
```bash
cd ~ && nf-core pipelines create --template-yaml ~/pneuma-template.yml
mv Julio92-C-pneuma pneuma 2>/dev/null || true
cd ~/pneuma && ls
```
Expected: `main.nf nextflow.config workflows subworkflows modules conf assets ...`. If the command rejects a `skip_features` key, delete that key and re-run.

- [ ] **Step 3: Verify the template runs**

Run: `nextflow run . -profile test,docker --outdir /tmp/pneuma_tpl`
Expected: `Pipeline completed successfully`.

- [ ] **Step 4: Add MIT license, spec, pytest harness**

`LICENSE`:
```
MIT License

Copyright (c) 2026 Julio Cesar Ortega Cambara

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
Copy the spec:
```bash
mkdir -p docs/superpowers/specs
cp /mnt/c/Users/julio/Desktop/Metagenomics_pipeline_automation/docs/superpowers/specs/2026-10-08-pneuma-design.md docs/superpowers/specs/
cp /mnt/c/Users/julio/Desktop/Metagenomics_pipeline_automation/docs/superpowers/plans/2026-10-08-pneuma-phase1-core.md docs/superpowers/
```
`tests/python/conftest.py`:
```python
import importlib.util
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
DATA = Path(__file__).resolve().parents[1] / "data"


def load_bin(name):
    """Import bin/<name>.py as a module (bin/ is not a package)."""
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module
```
`.github/workflows/pytest.yml`:
```yaml
name: pytest
on: [push, pull_request]
jobs:
  pytest:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: { python-version: "3.9" }
      - run: pip install pytest
      - run: pytest tests/python -v
```

- [ ] **Step 5: Create the private GitHub repo and push**

```bash
git add -A && git commit -m "chore: scaffold pneuma from nf-core template, add MIT license and design spec

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
gh repo create Julio92-C/pneuma --private --source . --push
gh repo view Julio92-C/pneuma --json visibility -q .visibility
```
Expected: `PRIVATE`.

---

### Task 2: Samplesheet schema, run grouping, control warning

**Files:**
- Modify: `assets/schema_input.json` (replace), `subworkflows/local/utils_nfcore_pneuma_pipeline/main.nf` (channel creation + helpers), `assets/samplesheet_test.csv` (replace)
- Test: `subworkflows/local/utils_nfcore_pneuma_pipeline/tests/main.nf.test`

**Interfaces:**
- Produces: `PIPELINE_INITIALISATION.out.samplesheet` emitting `[ meta, [fastq, ...] ]` with `meta = [id, batch, site, is_control (Boolean), run_ids (List), single_end: true]`. One element per sample; fastqs ordered by `run_id`.

- [ ] **Step 1: Write the schema**

`assets/schema_input.json`:
```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://raw.githubusercontent.com/Julio92-C/pneuma/main/assets/schema_input.json",
  "title": "PNEUMA samplesheet",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "sample": {
        "type": "string",
        "pattern": "^[A-Za-z0-9][A-Za-z0-9-]*$",
        "errorMessage": "sample must match ^[A-Za-z0-9][A-Za-z0-9-]*$ (letters, digits, '-'; no '_' or '.')",
        "meta": ["id"]
      },
      "run_id": { "type": "string", "pattern": "^\\S+$", "errorMessage": "run_id must be non-empty, no spaces", "meta": ["run_id"] },
      "batch":  { "type": "string", "pattern": "^\\S+$", "errorMessage": "batch must be non-empty, no spaces", "meta": ["batch"] },
      "site":   { "type": "string", "meta": ["site"] },
      "fastq": {
        "type": "string",
        "format": "file-path",
        "exists": true,
        "pattern": "^\\S+\\.f(ast)?q\\.gz$",
        "errorMessage": "fastq must be an existing .fastq.gz/.fq.gz path"
      },
      "is_control": { "type": "boolean", "meta": ["is_control"] }
    },
    "required": ["sample", "run_id", "batch", "fastq", "is_control"]
  },
  "uniqueEntries": ["sample", "run_id"]
}
```

- [ ] **Step 2: Write the test samplesheet**

`assets/samplesheet_test.csv`:
```
sample,run_id,batch,site,fastq,is_control
EBF0-0001,RUN1,TEST,SITEA,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz,false
EBF0-0001,RUN2,TEST,SITEA,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz,false
SPF0-0001,RUN1,TEST,SITEA,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz,false
CNF0-0001,RUN1,TEST,SITEA,https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz,true
```

- [ ] **Step 3: Write the failing nf-test**

`subworkflows/local/utils_nfcore_pneuma_pipeline/tests/main.nf.test`:
```groovy
nextflow_workflow {
    name "Test PIPELINE_INITIALISATION"
    script "../main.nf"
    workflow "PIPELINE_INITIALISATION"

    test("groups runs per sample and keeps control flag") {
        when {
            params { outdir = "$outputDir" }
            workflow {
                """
                input[0] = false
                input[1] = false
                input[2] = false
                input[3] = []
                input[4] = "$outputDir"
                input[5] = "${projectDir}/assets/samplesheet_test.csv"
                """
            }
        }
        then {
            assert workflow.success
            def items = workflow.out.samplesheet
            assert items.size() == 3
            def eb = items.find { it[0].id == 'EBF0-0001' }
            assert eb[0].run_ids == ['RUN1', 'RUN2']
            assert eb[1].size() == 2
            assert eb[0].single_end == true
            def cn = items.find { it[0].id == 'CNF0-0001' }
            assert cn[0].is_control == true
        }
    }
}
```
Open the template's `subworkflows/local/utils_nfcore_pneuma_pipeline/main.nf` and check the `take:` block of `PIPELINE_INITIALISATION`; reorder the six `input[n]` lines above so they match it (template order is `version, validate_params, monochrome_logs, nextflow_cli_args, outdir, input`).

- [ ] **Step 4: Run it to confirm it fails**

Run: `nf-test test subworkflows/local/utils_nfcore_pneuma_pipeline/tests/main.nf.test`
Expected: FAIL — template still expects `fastq_1/fastq_2`.

- [ ] **Step 5: Replace channel creation and add helpers**

In `subworkflows/local/utils_nfcore_pneuma_pipeline/main.nf`, replace the template's `Channel.fromList(samplesheetToList(...)) ... .set { ch_samplesheet }` block with:
```groovy
    def rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    warnRunsWithoutControls(rows)

    Channel
        .fromList(rows)
        .map { meta, fastq -> [ meta.id, meta, fastq ] }
        .groupTuple()
        .map { id, metas, fastqs -> groupRuns(id, metas, fastqs) }
        .set { ch_samplesheet }
```
Delete the template's `validateInputSamplesheet` function and add at the bottom of the file:
```groovy
//
// Merge all runs of one sample into a single meta + ordered fastq list.
// batch / site / is_control must agree across a sample's runs.
//
def groupRuns(id, metas, fastqs) {
    ['batch', 'site', 'is_control'].each { key ->
        def values = metas.collect { it[key] }.unique()
        if (values.size() > 1) {
            error("Sample '${id}' has conflicting '${key}' values across runs: ${values}")
        }
    }
    def order = (0..<metas.size()).sort { metas[it].run_id }
    def meta = [
        id         : id,
        batch      : metas[0].batch,
        site       : metas[0].site ?: '',
        is_control : metas[0].is_control as Boolean,
        run_ids    : order.collect { metas[it].run_id },
        single_end : true
    ]
    return [ meta, order.collect { fastqs[it] } ]
}

//
// Warn (not fail) for every run_id that has no negative control.
//
def warnRunsWithoutControls(rows) {
    def runs = rows.collect { it[0].run_id }.unique().sort()
    def withCtrl = rows.findAll { it[0].is_control }.collect { it[0].run_id }.unique()
    def missing = runs - withCtrl
    if (missing) {
        log.warn "No negative control declared for run(s): ${missing.join(', ')}. " +
                 "Contamination assessment for these runs relies on controls from other runs."
    }
}
```

- [ ] **Step 6: Run the test to confirm it passes**

Run: `nf-test test subworkflows/local/utils_nfcore_pneuma_pipeline/tests/main.nf.test`
Expected: PASS, and the log contains `No negative control declared for run(s): RUN2`.

- [ ] **Step 7: Commit**

```bash
git add assets/schema_input.json assets/samplesheet_test.csv subworkflows/local/utils_nfcore_pneuma_pipeline
git commit -m "feat: samplesheet schema with per-run rows grouped per sample

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: `tag_reads_taxid.py` and `merge_abricate.py`

**Files:**
- Create: `bin/tag_reads_taxid.py`, `bin/merge_abricate.py`, `tests/data/classified_small.fastq`, `tests/data/abricate_S1.csv`, `tests/data/abricate_S2_empty.csv`
- Test: `tests/python/test_tag_reads_taxid.py`, `tests/python/test_merge_abricate.py`

**Interfaces:**
- Produces: CLI `tag_reads_taxid.py <in.fastq[.gz]> <out.fasta>`; function `tag_records(lines) -> Iterator[str]`.
- Produces: CLI `merge_abricate.py --out <csv> <in.csv>...`; constant `HEADER`; function `merge(paths, out_handle) -> int` (data rows written).

- [ ] **Step 1: Write fixtures**

`tests/data/classified_small.fastq` (the 2nd record's quality line deliberately starts with `@`):
```
@read1 runid=abc kraken:taxid|1301
ACGTACGT
+
IIIIIIII
@read2 kraken:taxid|727
GGCC
+
@@II
```
`tests/data/abricate_S1.csv`:
```
#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE
S1_kraken2ID_PFP.fasta,read1_kraken:taxid|1301,1,100,+,patB,1-100/100,========,0/0,100.00,99.00,card,X1,PatB transporter,fluoroquinolone
#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE
S1_kraken2ID_PFP.fasta,read9_kraken:taxid|727,5,90,-,IncFII,1-86/86,========,0/0,100.00,98.00,plasmidfinder,X2,IncFII,
```
`tests/data/abricate_S2_empty.csv`: an empty file (0 bytes) — create with `: > tests/data/abricate_S2_empty.csv`.

- [ ] **Step 2: Write the failing tests**

`tests/python/test_tag_reads_taxid.py`:
```python
import pytest
from conftest import DATA, load_bin

tag = load_bin("tag_reads_taxid")


def test_tags_reads_with_taxid():
    with open(DATA / "classified_small.fastq") as fh:
        out = list(tag.tag_records(fh))
    assert out == [
        ">read1_kraken:taxid|1301\nACGTACGT\n",
        ">read2_kraken:taxid|727\nGGCC\n",
    ]


def test_missing_taxid_raises():
    with pytest.raises(ValueError, match="kraken:taxid"):
        list(tag.tag_records(["@r1 nothing\n", "AC\n", "+\n", "II\n"]))


def test_truncated_record_raises():
    with pytest.raises(ValueError, match="truncated"):
        list(tag.tag_records(["@r1 kraken:taxid|1\n", "AC\n"]))


def test_cli_writes_fasta(tmp_path):
    out = tmp_path / "S1_kraken2ID_PFP.fasta"
    tag.main([str(DATA / "classified_small.fastq"), str(out)])
    assert out.read_text().startswith(">read1_kraken:taxid|1301\n")
```
`tests/python/test_merge_abricate.py`:
```python
import io
from conftest import DATA, load_bin

merge = load_bin("merge_abricate")


def test_single_header_and_rows_kept():
    buf = io.StringIO()
    n = merge.merge([DATA / "abricate_S1.csv", DATA / "abricate_S2_empty.csv"], buf)
    lines = buf.getvalue().splitlines()
    assert n == 2
    assert lines[0] == merge.HEADER
    assert sum(1 for l in lines if l.startswith("#FILE")) == 1
    assert lines[2].endswith("IncFII,")


def test_only_empty_inputs_gives_header_only():
    buf = io.StringIO()
    assert merge.merge([DATA / "abricate_S2_empty.csv"], buf) == 0
    assert buf.getvalue() == merge.HEADER + "\n"


def test_unexpected_header_raises(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("#FILE,SEQUENCE,GENE\nx,y,z\n")
    try:
        merge.merge([bad], io.StringIO())
    except ValueError as e:
        assert "header" in str(e)
    else:
        raise AssertionError("expected ValueError")
```

- [ ] **Step 3: Run tests to confirm they fail**

Run: `pytest tests/python -v`
Expected: FAIL — `bin/tag_reads_taxid.py` not found.

- [ ] **Step 4: Implement**

`bin/tag_reads_taxid.py`:
```python
#!/usr/bin/env python3
"""Convert Kraken2 --classified-out FASTQ into taxid-tagged FASTA for ABRicate.

Input header : @<read_id> [...] kraken:taxid|<taxid> [...]
Output header: ><read_id>_kraken:taxid|<taxid>

This is the read-naming contract MERIDIAN uses to link every ABRicate hit
to the taxon of the read that carries it.
"""
import argparse
import gzip
import re
import sys

TAXID_RE = re.compile(r"kraken:taxid\|(\d+)")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def tag_records(lines):
    it = iter(lines)
    for header in it:
        header = header.rstrip("\r\n")
        if not header:
            continue
        try:
            seq = next(it).rstrip("\r\n")
            next(it)
            next(it)
        except StopIteration:
            raise ValueError(f"truncated FASTQ record after header: {header[:60]}")
        if not header.startswith("@"):
            raise ValueError(f"not a FASTQ header: {header[:60]}")
        match = TAXID_RE.search(header)
        if match is None:
            raise ValueError(f"no kraken:taxid tag in header: {header[:60]}")
        read_id = header[1:].split()[0]
        yield f">{read_id}_kraken:taxid|{match.group(1)}\n{seq}\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq", help="Kraken2 classified reads (.fastq or .fastq.gz)")
    parser.add_argument("output", help="output FASTA")
    args = parser.parse_args(argv)
    n = 0
    with open_text(args.fastq) as fh, open(args.output, "w") as out:
        for record in tag_records(fh):
            out.write(record)
            n += 1
    print(f"tagged {n} reads", file=sys.stderr)


if __name__ == "__main__":
    main()
```
`bin/merge_abricate.py`:
```python
#!/usr/bin/env python3
"""Concatenate per-sample ABRicate --csv outputs into one summary_reportClean.csv.

Lines are copied verbatim (ABRicate does not quote fields, so re-parsing
with csv could split PRODUCT values). Repeated headers and blank lines are
dropped; one header is written first.
"""
import argparse
import sys

HEADER = ("#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,"
          "%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE")


def merge(paths, out):
    out.write(HEADER + "\n")
    rows = 0
    for path in sorted(str(p) for p in paths):
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#FILE"):
                    if line != HEADER:
                        raise ValueError(f"unexpected ABRicate header in {path}: {line}")
                    continue
                out.write(line + "\n")
                rows += 1
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args(argv)
    with open(args.out, "w") as out:
        n = merge(args.inputs, out)
    print(f"merged {n} ABRicate hits from {len(args.inputs)} files", file=sys.stderr)


if __name__ == "__main__":
    main()
```
Then: `chmod +x bin/tag_reads_taxid.py bin/merge_abricate.py`

- [ ] **Step 5: Run tests to confirm they pass**

Run: `pytest tests/python -v`
Expected: 7 passed.

- [ ] **Step 6: Commit**

```bash
git add bin/tag_reads_taxid.py bin/merge_abricate.py tests/python tests/data
git commit -m "feat: taxid read tagging and ABRicate merge helpers

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: `arrange_reports.py` (Kraken2/Bracken tables for MERIDIAN)

**Files:**
- Create: `bin/arrange_reports.py`, `tests/data/bracken_combined.txt`, `tests/data/kraken_combined.txt`
- Test: `tests/python/test_arrange_reports.py`

**Interfaces:**
- Produces: CLI `arrange_reports.py --kraken <combined.txt> --bracken <combined.txt> --kraken-out kraken2_combined_reports.txt --bracken-out bracken_arranged.csv`; functions `read_from_perc(lines) -> list[str]`, `arrange_bracken(lines) -> list[list[str]]` (header row first).
- Behaviour mirrors `projects/lung_microbiome_2026/prep_inputs.R` §1–2 (MERIDIAN repo): strip the `combine_kreports.py` preamble; Bracken wide table `#perc,tot_all,<samples...>,taxid,name` where `#perc = round(100*tot_all / sum(tot_all of root+unclassified), 4)`.

- [ ] **Step 1: Write fixtures** (tab-separated; `\t` below means a literal TAB)

`tests/data/bracken_combined.txt`:
```
#Number of Samples: 2
#Total Number of Reads: 100
#Sample 1: S1
#Sample 2: S2
#perc\ttot_all\ttot_lvl\tS1_all\tS1_lvl\tS2_all\tS2_lvl\tlvl_type\ttaxid\tname
10.00\t10\t10\t4\t4\t6\t6\tU\t0\tunclassified
90.00\t90\t0\t36\t0\t54\t0\tR\t1\troot
60.00\t60\t60\t30\t30\t30\t30\tS\t1301\t    Streptococcus pneumoniae
30.00\t30\t30\t6\t6\t24\t24\tS\t727\t    Haemophilus influenzae
```
`tests/data/kraken_combined.txt`: the same content as above (only the preamble stripping is tested on it).

Create both with `printf` so tabs are real, e.g.:
```bash
printf '#Number of Samples: 2\n#Total Number of Reads: 100\n#Sample 1: S1\n#Sample 2: S2\n#perc\ttot_all\ttot_lvl\tS1_all\tS1_lvl\tS2_all\tS2_lvl\tlvl_type\ttaxid\tname\n10.00\t10\t10\t4\t4\t6\t6\tU\t0\tunclassified\n90.00\t90\t0\t36\t0\t54\t0\tR\t1\troot\n60.00\t60\t60\t30\t30\t30\t30\tS\t1301\t    Streptococcus pneumoniae\n30.00\t30\t30\t6\t6\t24\t24\tS\t727\t    Haemophilus influenzae\n' > tests/data/bracken_combined.txt
cp tests/data/bracken_combined.txt tests/data/kraken_combined.txt
```

- [ ] **Step 2: Write the failing tests**

`tests/python/test_arrange_reports.py`:
```python
import csv
import pytest
from conftest import DATA, load_bin

ar = load_bin("arrange_reports")


def lines(name):
    return (DATA / name).read_text().splitlines()


def test_read_from_perc_strips_preamble():
    out = ar.read_from_perc(lines("kraken_combined.txt"))
    assert out[0].startswith("#perc\t")
    assert len(out) == 5


def test_read_from_perc_requires_one_header():
    with pytest.raises(ValueError):
        ar.read_from_perc(["#Number of Samples: 1", "no header here"])


def test_arrange_bracken_shape_and_values():
    rows = ar.arrange_bracken(ar.read_from_perc(lines("bracken_combined.txt")))
    assert rows[0] == ["#perc", "tot_all", "S1", "S2", "taxid", "name"]
    by_name = {r[-1]: r for r in rows[1:]}
    assert by_name["Streptococcus pneumoniae"][1:4] == ["60", "30", "30"]
    # denominator = root(90) + unclassified(10) = 100
    assert by_name["Streptococcus pneumoniae"][0] == "60"
    assert by_name["Haemophilus influenzae"][0] == "30"
    assert by_name["Haemophilus influenzae"][4] == "727"


def test_cli_writes_both_files(tmp_path):
    k_out, b_out = tmp_path / "k.txt", tmp_path / "b.csv"
    ar.main(["--kraken", str(DATA / "kraken_combined.txt"),
             "--bracken", str(DATA / "bracken_combined.txt"),
             "--kraken-out", str(k_out), "--bracken-out", str(b_out)])
    assert k_out.read_text().startswith("#perc\t")
    with open(b_out) as fh:
        assert next(csv.reader(fh))[0] == "#perc"
```

- [ ] **Step 3: Run tests to confirm they fail**

Run: `pytest tests/python/test_arrange_reports.py -v`
Expected: FAIL — module not found.

- [ ] **Step 4: Implement**

`bin/arrange_reports.py`:
```python
#!/usr/bin/env python3
"""Turn KrakenTools combine_kreports.py outputs into MERIDIAN inputs.

--kraken  combined Kraken2 kreport  -> preamble stripped (starts at '#perc')
--bracken combined Bracken kreport  -> bracken_arranged.csv:
          #perc,tot_all,<sample...>,taxid,name
          #perc = 100 * tot_all / (tot_all[root] + tot_all[unclassified]), 4 dp
"""
import argparse
import csv


def read_from_perc(lines):
    starts = [i for i, line in enumerate(lines) if line.startswith("#perc\t")]
    if len(starts) != 1:
        raise ValueError(f"expected exactly one '#perc' header line, found {len(starts)}")
    return [line.rstrip("\r\n") for line in lines[starts[0]:]]


def fmt_number(value):
    value = round(value, 4)
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.4f}".rstrip("0").rstrip(".")


def arrange_bracken(lines):
    header = lines[0].split("\t")
    sample_cols = [c for c in header if c.endswith("_all") and c != "tot_all"]
    samples = [c[: -len("_all")] for c in sample_cols]
    idx = {c: header.index(c) for c in sample_cols + ["taxid", "name"]}
    body = []
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        counts = [float(fields[idx[c]]) for c in sample_cols]
        body.append((counts, fields[idx["taxid"]], fields[idx["name"]].strip()))
    n_total = sum(sum(c) for c, _, name in body if name in ("root", "unclassified"))
    if n_total == 0:
        raise ValueError("root + unclassified total is zero; is this a Bracken kreport?")
    rows = [["#perc", "tot_all"] + samples + ["taxid", "name"]]
    for counts, taxid, name in body:
        tot = sum(counts)
        rows.append([fmt_number(100 * tot / n_total), fmt_number(tot)]
                    + [fmt_number(c) for c in counts] + [taxid, name])
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kraken", required=True)
    parser.add_argument("--bracken", required=True)
    parser.add_argument("--kraken-out", required=True)
    parser.add_argument("--bracken-out", required=True)
    args = parser.parse_args(argv)

    with open(args.kraken) as fh:
        kraken = read_from_perc(fh.readlines())
    with open(args.kraken_out, "w") as out:
        out.write("\n".join(kraken) + "\n")

    with open(args.bracken) as fh:
        rows = arrange_bracken(read_from_perc(fh.readlines()))
    with open(args.bracken_out, "w", newline="") as out:
        csv.writer(out, lineterminator="\n").writerows(rows)


if __name__ == "__main__":
    main()
```
Then `chmod +x bin/arrange_reports.py`.

- [ ] **Step 5: Run tests to confirm they pass**

Run: `pytest tests/python -v`
Expected: all pass (11).

- [ ] **Step 6: Commit**

```bash
git add bin/arrange_reports.py tests/python/test_arrange_reports.py tests/data/*_combined.txt
git commit -m "feat: arrange combined Kraken2/Bracken reports into MERIDIAN inputs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: `make_meridian_config.py` + MERIDIAN config template

**Files:**
- Create: `bin/make_meridian_config.py`, `assets/meridian_config.template.yaml`, `tests/data/metadata_test.csv`
- Test: `tests/python/test_make_meridian_config.py`

**Interfaces:**
- Produces: CLI `make_meridian_config.py --template T --out config.yaml --study-id ID --metadata meta.csv --group COL [--random-effect COL] [--controls A,B]`; function `render(template_text, study_id, metadata_name, group, random_effect, controls) -> str`; function `check_metadata(path, group, random_effect, controls) -> None` (raises `ValueError`).
- The rendered config keeps the literal token `__PROJECT_ROOT__`; the MERIDIAN process replaces it with its absolute work dir (Task 10).

- [ ] **Step 1: Write the template**

`assets/meridian_config.template.yaml` (derived from MERIDIAN v1.0.0 `config/config_template.yaml`; `alpha_diversity` added, `diversity_csv` removed):
```yaml
# Generated by PNEUMA bin/make_meridian_config.py — do not edit by hand.
study:
  id: {{STUDY_ID}}
  name: "{{STUDY_ID}} (PNEUMA)"
  description: >
    MERIDIAN cohort analysis generated by the PNEUMA pipeline.
project_root: "__PROJECT_ROOT__"
inputs:
  abricate_summary: "inputs/summary_reportClean.csv"
  kraken2_report:   "inputs/kraken2_combined_reports.txt"
  bracken_report:   "inputs/bracken_arranged.csv"
  recentrifuge_csv: "inputs/pneuma.rcf.data.csv"
metadata:
  file: "inputs/{{METADATA_FILE}}"
  sample_id_col: "sample"
  group_cols:
    - {{GROUP_COL}}
  random_effect: {{RANDOM_EFFECT}}
  fixed_effects:
    - {{GROUP_COL}}
  controls: {{CONTROLS}}
taxid_fixes_file: ~
alpha_diversity:
  source: computed
filters:
  min_count_per_sample:   17
  min_count_for_species:  500
  min_prevalence_aldex:   2
stages:
  clean_data:             true
  normalisation:          true
  taxonomy:               true
  relative_abundance:     true
  alpha_diversity:        true
  beta_diversity:         true
  differential_abundance: true
  resistome:              true
  virulome:               true
  mobilome:               true
  network:                true
  report:                 true
  panels:                 true
stats:
  permanova_permutations: 9999
  aldex_mc_samples:       128
  alpha: 0.05
  padjust_method: BH
report:
  top_n_taxa:              25
  rank_taxa_by:            prevalence_x_abundance
  network_top_nodes:       50
  network_edge_min_weight: 0.3
  show_per_group_chord:    true
  show_per_group_sankey:   true
network:
  sankey:
    top_genes:  30
    rank_by:    prevalence_x_abundance
panels:
  composite_unified:     true
  composite_base_size:   8
  composite_legend_size: 8
  composite_margin_pt:   6
  strip_plot_title:      true
  show_title:            false
  png:                   true
  tiff:                  false
  pdf:                   false
outputs:
  datasets_dir: "Datasets"
  figures_dir:  "Figures"
  report_html:  "Reports/pipeline_report.html"
  log_file:     "Reports/pipeline.log"
```
Before committing, diff the non-placeholder keys against MERIDIAN v1.0.0:
```bash
git -C /mnt/c/Users/julio/Desktop/Metagenomics_pipeline_automation show v1.0.0:config/config_template.yaml | grep -vE '^\s*#|^\s*$' > /tmp/m.yaml
diff <(grep -oE '^\s*[a-z_0-9]+:' /tmp/m.yaml) <(grep -oE '^\s*[a-z_0-9]+:' assets/meridian_config.template.yaml)
```
Expected differences only: `diversity_csv` removed, `alpha_diversity`/`source` added.

`tests/data/metadata_test.csv`:
```
sample,sample_type,ID
EBF0-0001,EB,0001
SPF0-0001,SP,0001
CNF0-0001,CONTROL,0001
```

- [ ] **Step 2: Write the failing tests**

`tests/python/test_make_meridian_config.py`:
```python
import pytest
from conftest import BIN, DATA, load_bin

mm = load_bin("make_meridian_config")
TEMPLATE = (BIN.parent / "assets" / "meridian_config.template.yaml").read_text()


def test_render_fills_all_placeholders():
    text = mm.render(TEMPLATE, "lung", "metadata_test.csv", "sample_type", "ID", ["CNF0-0001"])
    assert "{{" not in text
    assert "group_cols:\n    - sample_type" in text
    assert "random_effect: ID" in text
    assert 'controls: ["CNF0-0001"]' in text
    assert 'file: "inputs/metadata_test.csv"' in text
    assert "__PROJECT_ROOT__" in text


def test_render_defaults():
    text = mm.render(TEMPLATE, "lung", "m.csv", "sample_type", None, [])
    assert "random_effect: ~" in text
    assert "controls: []" in text


def test_check_metadata_ok():
    mm.check_metadata(DATA / "metadata_test.csv", "sample_type", "ID", ["CNF0-0001"])


@pytest.mark.parametrize("group,re_col,controls,msg", [
    ("missing_col", None, [], "group column"),
    ("sample_type", "nope", [], "random-effect column"),
    ("sample_type", None, ["NOT-THERE"], "control"),
])
def test_check_metadata_errors(group, re_col, controls, msg):
    with pytest.raises(ValueError, match=msg):
        mm.check_metadata(DATA / "metadata_test.csv", group, re_col, controls)


def test_cli(tmp_path):
    out = tmp_path / "config.yaml"
    mm.main(["--template", str(BIN.parent / "assets" / "meridian_config.template.yaml"),
             "--out", str(out), "--study-id", "lung",
             "--metadata", str(DATA / "metadata_test.csv"),
             "--group", "sample_type", "--controls", "CNF0-0001"])
    assert "controls: [\"CNF0-0001\"]" in out.read_text()
```

- [ ] **Step 3: Run tests to confirm they fail**

Run: `pytest tests/python/test_make_meridian_config.py -v`
Expected: FAIL — module not found.

- [ ] **Step 4: Implement**

`bin/make_meridian_config.py`:
```python
#!/usr/bin/env python3
"""Render a MERIDIAN v1.0.0 config.yaml from PNEUMA's template.

Placeholders: {{STUDY_ID}} {{METADATA_FILE}} {{GROUP_COL}} {{RANDOM_EFFECT}} {{CONTROLS}}.
'__PROJECT_ROOT__' is left for the MERIDIAN process to replace with its work dir.
"""
import argparse
import csv
from pathlib import Path


def yaml_list(items):
    return "[" + ", ".join(f'"{i}"' for i in items) + "]"


def render(template_text, study_id, metadata_name, group, random_effect, controls):
    text = (template_text
            .replace("{{STUDY_ID}}", study_id)
            .replace("{{METADATA_FILE}}", metadata_name)
            .replace("{{GROUP_COL}}", group)
            .replace("{{RANDOM_EFFECT}}", random_effect or "~")
            .replace("{{CONTROLS}}", yaml_list(controls)))
    if "{{" in text:
        raise ValueError("unfilled placeholder left in MERIDIAN config template")
    return text


def check_metadata(path, group, random_effect, controls):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames or []
        samples = {row["sample"] for row in reader} if "sample" in columns else set()
    if "sample" not in columns:
        raise ValueError(f"metadata {path} has no 'sample' column")
    if group not in columns:
        raise ValueError(f"group column '{group}' not in metadata columns {columns}")
    if random_effect and random_effect not in columns:
        raise ValueError(f"random-effect column '{random_effect}' not in metadata columns {columns}")
    missing = [c for c in controls if c not in samples]
    if missing:
        raise ValueError(f"control sample(s) not in metadata: {missing}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--random-effect", default=None)
    parser.add_argument("--controls", default="")
    args = parser.parse_args(argv)

    controls = [c for c in args.controls.split(",") if c]
    check_metadata(args.metadata, args.group, args.random_effect, controls)
    text = render(Path(args.template).read_text(), args.study_id,
                  Path(args.metadata).name, args.group, args.random_effect, controls)
    Path(args.out).write_text(text)


if __name__ == "__main__":
    main()
```
Then `chmod +x bin/make_meridian_config.py`.

- [ ] **Step 5: Run tests to confirm they pass**

Run: `pytest tests/python -v`
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add bin/make_meridian_config.py assets/meridian_config.template.yaml tests/python/test_make_meridian_config.py tests/data/metadata_test.csv
git commit -m "feat: MERIDIAN config rendering with metadata checks

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: PREPROCESS subworkflow

**Files:**
- Create: `modules/local/hostile_clean/main.nf`, `modules/local/hostile_clean/tests/main.nf.test`, `subworkflows/local/preprocess/main.nf`, `subworkflows/local/preprocess/tests/main.nf.test`
- Install: nf-core `cat/fastq`, `porechop/abi`, `nanoq`, `nanoplot`, `seqkit/stats`
- Modify: `nextflow.config` (params), `conf/modules.config`, `conf/test.config`

**Interfaces:**
- Consumes: `[ meta, [fastq...] ]` from Task 2.
- Produces: `PREPROCESS.out.reads` → `[ meta, clean.fastq.gz ]` with `meta.low_depth` (Boolean) and `meta.n_reads` (Long) added; `PREPROCESS.out.multiqc` (files); `PREPROCESS.out.versions`.

- [ ] **Step 1: Install nf-core modules and read their interfaces**

```bash
for m in cat/fastq porechop/abi nanoq nanoplot seqkit/stats; do nf-core modules install $m; done
grep -A6 "input:" modules/nf-core/{cat/fastq,porechop/abi,nanoq,nanoplot,seqkit/stats}/main.nf
```
Expected signatures (adjust the subworkflow calls in Step 4 if any differ):
- `CAT_FASTQ(tuple val(meta), path(reads, stageAs: "input*/*"))` → `reads`
- `PORECHOP_ABI(tuple val(meta), path(reads), path(custom_adapters))` → `reads` (pass `[]` for custom_adapters; if the installed version has no second input, drop it)
- `NANOQ(tuple val(meta), path(ontreads), val(output_format))` → `reads`, `stats`
- `NANOPLOT(tuple val(meta), path(ontfile))` → `txt`, `html`
- `SEQKIT_STATS(tuple val(meta), path(reads))` → `stats`

- [ ] **Step 2: Pin the Hostile container and confirm its CLI**

```bash
curl -s "https://quay.io/api/v1/repository/biocontainers/hostile/tag/?limit=5&onlyActiveTags=true" | python3 -c "import sys,json;[print(t['name']) for t in json.load(sys.stdin)['tags']]"
```
Pick the newest tag `T` (≥ 2.0). Then:
```bash
docker run --rm quay.io/biocontainers/hostile:T hostile clean --help | grep -E -- "--fastq1|--aligner|--index|--airplane|--output|--threads"
```
Expected: all six flags listed. If `--airplane` is absent, delete it from the script in Step 3.

- [ ] **Step 3: Write HOSTILE_CLEAN and its test**

`modules/local/hostile_clean/main.nf` (replace `T` with the tag from Step 2):
```groovy
process HOSTILE_CLEAN {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/hostile:T'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.hostile.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.hostile.json")    , emit: log
    path "versions.yml"                                  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    hostile clean \\
        --fastq1 ${reads} \\
        --aligner minimap2 \\
        --index ${index} \\
        --threads ${task.cpus} \\
        --output hostile_out \\
        --airplane \\
        ${args} > ${meta.id}.hostile.json
    mv hostile_out/*.fastq.gz ${meta.id}.hostile.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hostile: \$(hostile --version 2>&1 | sed 's/^.*hostile //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    echo -n | gzip > ${meta.id}.hostile.fastq.gz
    echo '{}' > ${meta.id}.hostile.json
    echo '"${task.process}": {hostile: stub}' > versions.yml
    """
}
```
`modules/local/hostile_clean/tests/main.nf.test`:
```groovy
nextflow_process {
    name "Test HOSTILE_CLEAN"
    script "../main.nf"
    process "HOSTILE_CLEAN"

    test("removes reads against a small host fasta") {
        when {
            process {
                """
                input[0] = [ [id:'S1'], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz', checkIfExists: true) ]
                input[1] = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta', checkIfExists: true)
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.reads[0][1]).linesGzip.size() > 0
            assert path(process.out.log[0][1]).text.contains('reads_in')
        }
    }
}
```
Run: `nf-test test modules/local/hostile_clean/tests/main.nf.test --profile docker`
Expected: PASS. (SARS-CoV-2 reads do not map to the human slice, so reads survive.) If the JSON has no `reads_in` key, open it, pick the input-read-count key it does have, and update the assertion.

- [ ] **Step 4: Write the PREPROCESS subworkflow**

`subworkflows/local/preprocess/main.nf`:
```groovy
include { CAT_FASTQ     } from '../../../modules/nf-core/cat/fastq/main'
include { PORECHOP_ABI  } from '../../../modules/nf-core/porechop/abi/main'
include { NANOQ         } from '../../../modules/nf-core/nanoq/main'
include { HOSTILE_CLEAN } from '../../../modules/local/hostile_clean/main'
include { NANOPLOT      } from '../../../modules/nf-core/nanoplot/main'
include { SEQKIT_STATS  } from '../../../modules/nf-core/seqkit/stats/main'

workflow PREPROCESS {
    take:
    ch_reads      // [ meta, [fastq...] ]
    host_index    // path: fasta or .mmi

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    ch_reads
        .branch { meta, fastqs ->
            multi : fastqs.size() > 1
            single: true
                return [ meta, fastqs[0] ]
        }
        .set { ch_split }

    CAT_FASTQ(ch_split.multi)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())
    ch_merged   = CAT_FASTQ.out.reads.mix(ch_split.single)

    if (params.run_porechop) {
        PORECHOP_ABI(ch_merged, [])
        ch_merged   = PORECHOP_ABI.out.reads
        ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions.first())
    }

    NANOQ(ch_merged, 'fastq.gz')
    ch_versions = ch_versions.mix(NANOQ.out.versions.first())

    HOSTILE_CLEAN(NANOQ.out.reads, host_index)
    ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions.first())

    NANOPLOT(HOSTILE_CLEAN.out.reads)
    SEQKIT_STATS(HOSTILE_CLEAN.out.reads)
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first(), SEQKIT_STATS.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(NANOPLOT.out.txt.collect { it[1] }, NANOQ.out.stats.collect { it[1] })

    // seqkit stats -T: header + one row; column 3 (0-based) = num_seqs
    ch_counts = SEQKIT_STATS.out.stats.map { meta, tsv ->
        [ meta.id, tsv.readLines()[1].split('\t')[3].replace(',', '') as Long ]
    }

    HOSTILE_CLEAN.out.reads
        .map { meta, fq -> [ meta.id, meta, fq ] }
        .join(ch_counts)
        .map { id, meta, fq, n ->
            if (n < params.min_reads) {
                log.warn "Sample ${id}: ${n} reads after host removal (< --min_reads ${params.min_reads}); flagged LOW_DEPTH"
            }
            [ meta + [ n_reads: n, low_depth: n < params.min_reads ], fq ]
        }
        .set { ch_clean }

    emit:
    reads    = ch_clean
    multiqc  = ch_multiqc
    versions = ch_versions
}
```

- [ ] **Step 5: Add params and module options**

In `nextflow.config` `params { }` add:
```groovy
    // PREPROCESS
    host_index           = null
    run_porechop         = false
    nanoq_min_length     = 100     // initial value; revisited in Phase 3 regression
    nanoq_min_qscore     = 7       // initial value; revisited in Phase 3 regression
    min_reads            = 100
```
In `conf/modules.config` `process { }` add:
```groovy
    withName: 'NANOQ' {
        ext.args   = { "--min-len ${params.nanoq_min_length} --min-qual ${params.nanoq_min_qscore}" }
        ext.prefix = { "${meta.id}.nanoq" }
        publishDir = [ path: { "${params.outdir}/preprocess/nanoq" }, mode: params.publish_dir_mode, pattern: '*.stats' ]
    }
    withName: 'HOSTILE_CLEAN' {
        publishDir = [ path: { "${params.outdir}/preprocess/clean_reads" }, mode: params.publish_dir_mode ]
    }
    withName: 'NANOPLOT' {
        publishDir = [ path: { "${params.outdir}/preprocess/nanoplot/${meta.id}" }, mode: params.publish_dir_mode ]
    }
    withName: 'SEQKIT_STATS' {
        ext.args   = '-T -a'
        publishDir = [ path: { "${params.outdir}/preprocess/seqkit" }, mode: params.publish_dir_mode ]
    }
```
In `conf/test.config` `params { }` set:
```groovy
    input      = "${projectDir}/assets/samplesheet_test.csv"
    host_index = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta'
    min_reads  = 10
```

- [ ] **Step 6: Write and run the subworkflow test**

`subworkflows/local/preprocess/tests/main.nf.test`:
```groovy
nextflow_workflow {
    name "Test PREPROCESS"
    script "../main.nf"
    workflow "PREPROCESS"
    config "../../../../conf/test.config"

    test("merges runs, filters, removes host, flags depth") {
        when {
            params { outdir = "$outputDir"; min_reads = 10 }
            workflow {
                """
                def fq = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz', checkIfExists: true)
                input[0] = Channel.of(
                    [ [id:'EBF0-0001', single_end:true, is_control:false, run_ids:['RUN1','RUN2']], [fq, fq] ],
                    [ [id:'SPF0-0001', single_end:true, is_control:false, run_ids:['RUN1']], [fq] ]
                )
                input[1] = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta', checkIfExists: true)
                """
            }
        }
        then {
            assert workflow.success
            def out = workflow.out.reads.collectEntries { [ it[0].id, it[0] ] }
            assert out['EBF0-0001'].n_reads == 2 * out['SPF0-0001'].n_reads
            assert out['SPF0-0001'].low_depth == false
        }
    }
}
```
Run: `nf-test test subworkflows/local/preprocess/tests/main.nf.test --profile docker`
Expected: PASS (EBF0-0001 has exactly twice the reads because both "runs" are the same file).

- [ ] **Step 7: Commit**

```bash
git add modules.json modules/nf-core modules/local/hostile_clean subworkflows/local/preprocess nextflow.config conf/modules.config conf/test.config
git commit -m "feat: PREPROCESS subworkflow (run merge, nanoq, hostile, depth flag)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: TAXONOMY subworkflow

**Files:**
- Create: `modules/local/combine_kreports/main.nf`, `modules/local/arrange_reports/main.nf`, `modules/local/arrange_reports/tests/main.nf.test`, `subworkflows/local/taxonomy/main.nf`, `subworkflows/local/taxonomy/tests/main.nf.test`
- Install: nf-core `untar`, `kraken2/kraken2`, `bracken/bracken`, `krakentools/kreport2krona`, `krona/ktimporttext`
- Modify: `nextflow.config`, `conf/modules.config`, `conf/test.config`

**Interfaces:**
- Consumes: `PREPROCESS.out.reads` `[ meta(+low_depth), fastq.gz ]`.
- Produces:
  - `TAXONOMY.out.classified` → `[ meta, *.classified.fastq.gz ]` (all samples)
  - `TAXONOMY.out.assignments` → `[ meta, *.classifiedreads.txt ]` (per-read Kraken2 output, all samples)
  - `TAXONOMY.out.kraken_combined` → `path kraken2_combined_reports.txt`
  - `TAXONOMY.out.bracken_arranged` → `path bracken_arranged.csv`
  - `TAXONOMY.out.versions`

- [ ] **Step 1: Install modules and read interfaces**

```bash
for m in untar kraken2/kraken2 bracken/bracken krakentools/kreport2krona krona/ktimporttext; do nf-core modules install $m; done
grep -A8 "input:" modules/nf-core/{kraken2/kraken2,bracken/bracken}/main.nf
grep -A12 "output:" modules/nf-core/{kraken2/kraken2,bracken/bracken}/main.nf
```
Expected: `KRAKEN2_KRAKEN2(tuple val(meta), path(reads), path db, val save_output_fastqs, val save_reads_assignment)` with outputs `classified_reads_fastq`, `classified_reads_assignment`, `report`; `BRACKEN_BRACKEN(tuple val(meta), path(kraken_report), path database)` with outputs `reports` (`.tsv`) and `txt` (Bracken-style kreport). Use those emit names in Step 4; if the installed Bracken names its kreport emit differently, use the emit whose pattern ends in `report.txt`.

- [ ] **Step 2: Check the test Kraken2/Bracken DB**

```bash
curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/db/kraken2_bracken.tar.gz | tar -tzf - | grep kmer_distrib
```
Expected: one or more `database<N>mers.kmer_distrib`. Note the smallest `N`; it is the `bracken_read_length` for `conf/test.config`.

- [ ] **Step 3: Write COMBINE_KREPORTS and ARRANGE_REPORTS**

Confirm the KrakenTools flag first:
```bash
docker run --rm quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0 combine_kreports.py --help | grep -- "--sample-names"
```
Expected: the flag is listed. If the tag does not exist, take the newest from `https://quay.io/api/v1/repository/biocontainers/krakentools/tag/` and use it below.

`modules/local/combine_kreports/main.nf`:
```groovy
process COMBINE_KREPORTS {
    tag "$prefix"
    label 'process_single'
    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'

    input:
    val  prefix
    val  names      // List<String>, same order as reports
    path reports, stageAs: 'reports/*'

    output:
    path "${prefix}.txt", emit: txt
    path "versions.yml" , emit: versions

    script:
    """
    combine_kreports.py \\
        -r ${reports.join(' ')} \\
        -o ${prefix}.txt \\
        --sample-names ${names.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krakentools: 1.2
    END_VERSIONS
    """
}
```
`modules/local/arrange_reports/main.nf`:
```groovy
process ARRANGE_REPORTS {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path kraken_combined
    path bracken_combined

    output:
    path "kraken2_combined_reports.txt", emit: kraken
    path "bracken_arranged.csv"        , emit: bracken
    path "versions.yml"                , emit: versions

    script:
    """
    arrange_reports.py \\
        --kraken ${kraken_combined} \\
        --bracken ${bracken_combined} \\
        --kraken-out kraken2_combined_reports.txt \\
        --bracken-out bracken_arranged.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
```
`modules/local/arrange_reports/tests/main.nf.test`:
```groovy
nextflow_process {
    name "Test ARRANGE_REPORTS"
    script "../main.nf"
    process "ARRANGE_REPORTS"

    test("writes MERIDIAN tables") {
        when {
            process {
                """
                input[0] = file("${projectDir}/tests/data/kraken_combined.txt", checkIfExists: true)
                input[1] = file("${projectDir}/tests/data/bracken_combined.txt", checkIfExists: true)
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.bracken[0]).readLines()[0] == '#perc,tot_all,S1,S2,taxid,name'
            assert path(process.out.kraken[0]).readLines()[0].startsWith('#perc')
        }
    }
}
```
Run: `nf-test test modules/local/arrange_reports/tests/main.nf.test --profile docker`
Expected: PASS.

- [ ] **Step 4: Write the TAXONOMY subworkflow**

`subworkflows/local/taxonomy/main.nf`:
```groovy
include { UNTAR                    } from '../../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2          } from '../../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN          } from '../../../modules/nf-core/bracken/bracken/main'
include { KRAKENTOOLS_KREPORT2KRONA} from '../../../modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_KTIMPORTTEXT       } from '../../../modules/nf-core/krona/ktimporttext/main'
include { COMBINE_KREPORTS as COMBINE_KRAKEN  } from '../../../modules/local/combine_kreports/main'
include { COMBINE_KREPORTS as COMBINE_BRACKEN } from '../../../modules/local/combine_kreports/main'
include { ARRANGE_REPORTS          } from '../../../modules/local/arrange_reports/main'

workflow TAXONOMY {
    take:
    ch_reads     // [ meta, fastq.gz ]
    kraken2_db   // path: directory or .tar.gz

    main:
    ch_versions = Channel.empty()

    if (kraken2_db.toString().endsWith('.tar.gz')) {
        UNTAR(Channel.value([ [id:'kraken2_db'], file(kraken2_db, checkIfExists: true) ]))
        ch_db = UNTAR.out.untar.map { it[1] }.first()
    } else {
        ch_db = Channel.value(file(kraken2_db, checkIfExists: true))
    }

    KRAKEN2_KRAKEN2(ch_reads, ch_db, true, true)
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions.first())

    // Bracken is skipped for LOW_DEPTH samples (no reliable re-estimation)
    BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report.filter { meta, r -> !meta.low_depth }, ch_db)
    ch_versions = ch_versions.mix(BRACKEN_BRACKEN.out.versions.first())

    def collectSorted = { ch ->
        ch.toSortedList { a, b -> a[0].id <=> b[0].id }
          .map { items -> [ items.collect { it[0].id }, items.collect { it[1] } ] }
    }
    def k = collectSorted(KRAKEN2_KRAKEN2.out.report).multiMap { names, files -> names: names; files: files }
    def b = collectSorted(BRACKEN_BRACKEN.out.txt).multiMap { names, files -> names: names; files: files }

    COMBINE_KRAKEN('kraken2_combined', k.names, k.files)
    COMBINE_BRACKEN('bracken_combined', b.names, b.files)
    ch_versions = ch_versions.mix(COMBINE_KRAKEN.out.versions)

    ARRANGE_REPORTS(COMBINE_KRAKEN.out.txt, COMBINE_BRACKEN.out.txt)
    ch_versions = ch_versions.mix(ARRANGE_REPORTS.out.versions)

    KRAKENTOOLS_KREPORT2KRONA(KRAKEN2_KRAKEN2.out.report)
    KRONA_KTIMPORTTEXT(KRAKENTOOLS_KREPORT2KRONA.out.txt)
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
    classified       = KRAKEN2_KRAKEN2.out.classified_reads_fastq
    assignments      = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    reports          = KRAKEN2_KRAKEN2.out.report
    kraken_combined  = ARRANGE_REPORTS.out.kraken
    bracken_arranged = ARRANGE_REPORTS.out.bracken
    versions         = ch_versions
}
```

- [ ] **Step 5: Params and module config**

`nextflow.config` params:
```groovy
    // TAXONOMY
    kraken2_db           = null
    kraken2_confidence   = 0
    bracken_level        = 'S'
    bracken_threshold    = 10
    bracken_read_length  = 300
```
`conf/modules.config`:
```groovy
    withName: 'KRAKEN2_KRAKEN2' {
        ext.args   = { "--confidence ${params.kraken2_confidence}" }
        publishDir = [ path: { "${params.outdir}/taxonomy/kraken2" }, mode: params.publish_dir_mode,
                       pattern: '*.{report.txt,classifiedreads.txt}' ]
    }
    withName: 'BRACKEN_BRACKEN' {
        ext.args      = { "-l ${params.bracken_level} -t ${params.bracken_threshold} -r ${params.bracken_read_length}" }
        errorStrategy = { task.exitStatus == 1 ? 'ignore' : 'retry' }  // no taxa above threshold → sample has no Bracken estimate
        publishDir    = [ path: { "${params.outdir}/taxonomy/bracken" }, mode: params.publish_dir_mode ]
    }
    withName: 'ARRANGE_REPORTS|COMBINE_KRAKEN|COMBINE_BRACKEN' {
        publishDir = [ path: { "${params.outdir}/taxonomy" }, mode: params.publish_dir_mode ]
    }
    withName: 'KRONA_KTIMPORTTEXT' {
        publishDir = [ path: { "${params.outdir}/taxonomy/krona" }, mode: params.publish_dir_mode ]
    }
```
`conf/test.config` params (use the `N` from Step 2):
```groovy
    kraken2_db          = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/db/kraken2_bracken.tar.gz'
    bracken_read_length = N
```

- [ ] **Step 6: Write and run the subworkflow test**

`subworkflows/local/taxonomy/tests/main.nf.test`:
```groovy
nextflow_workflow {
    name "Test TAXONOMY"
    script "../main.nf"
    workflow "TAXONOMY"
    config "../../../../conf/test.config"

    test("classifies, re-estimates and arranges") {
        when {
            params { outdir = "$outputDir" }
            workflow {
                """
                def fq = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz', checkIfExists: true)
                input[0] = Channel.of(
                    [ [id:'EBF0-0001', single_end:true, low_depth:false], fq ],
                    [ [id:'SPF0-0001', single_end:true, low_depth:false], fq ]
                )
                input[1] = params.kraken2_db
                """
            }
        }
        then {
            assert workflow.success
            def header = path(workflow.out.bracken_arranged[0]).readLines()[0]
            assert header == '#perc,tot_all,EBF0-0001,SPF0-0001,taxid,name'
            assert workflow.out.assignments.size() == 2
        }
    }
}
```
Run: `nf-test test subworkflows/local/taxonomy/tests/main.nf.test --profile docker`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add modules.json modules/nf-core modules/local/{combine_kreports,arrange_reports} subworkflows/local/taxonomy nextflow.config conf
git commit -m "feat: TAXONOMY subworkflow (Kraken2, Bracken, combined tables, Krona)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: CONTAMINATION subworkflow (Recentrifuge)

**Files:**
- Create: `tests/make_test_taxdump.py`, `tests/data/taxdump/{nodes,names}.dmp` (generated), `modules/local/recentrifuge/main.nf`, `modules/local/recentrifuge/tests/main.nf.test`, `subworkflows/local/contamination/main.nf`
- Modify: `nextflow.config`, `conf/modules.config`, `conf/test.config`

**Interfaces:**
- Consumes: `TAXONOMY.out.assignments` `[ meta, <id>.kraken2.classifiedreads.txt ]`; taxdump directory.
- Produces: `CONTAMINATION.out.rcf_csv` → `path pneuma.rcf.data.csv`; `CONTAMINATION.out.html` (both HTMLs); `versions`.

- [ ] **Step 1: Generate a test taxdump from a test Kraken2 report**

`tests/make_test_taxdump.py`:
```python
#!/usr/bin/env python3
"""Build a minimal NCBI-style taxdump (nodes.dmp, names.dmp) from a Kraken2 kreport.

Test fixture only: the kreport's indentation gives each taxon's parent.
Usage: make_test_taxdump.py <kreport.txt> <outdir>
"""
import sys
from pathlib import Path

RANKS = {"R": "no rank", "D": "superkingdom", "K": "kingdom", "P": "phylum", "C": "class",
         "O": "order", "F": "family", "G": "genus", "S": "species"}


def main(report, outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    stack = []  # (depth, taxid)
    nodes, names = [], []
    for line in Path(report).read_text().splitlines():
        f = line.split("\t")
        if len(f) < 6 or f[4] == "0":
            continue
        rank_code, taxid, raw = f[3], f[4], f[5]
        depth = (len(raw) - len(raw.lstrip(" "))) // 2
        while stack and stack[-1][0] >= depth:
            stack.pop()
        parent = stack[-1][1] if stack else "1"
        rank = RANKS.get(rank_code[0], "no rank")
        nodes.append(f"{taxid}\t|\t{parent if taxid != '1' else '1'}\t|\t{rank}\t|")
        names.append(f"{taxid}\t|\t{raw.strip()}\t|\t\t|\tscientific name\t|")
        stack.append((depth, taxid))
    (out / "nodes.dmp").write_text("\n".join(nodes) + "\n")
    (out / "names.dmp").write_text("\n".join(names) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
```
Generate it from the report produced by the Task 7 test:
```bash
REPORT=$(find .nf-test -name 'EBF0-0001*.kraken2.report.txt' | head -1)
python3 tests/make_test_taxdump.py "$REPORT" tests/data/taxdump
head -3 tests/data/taxdump/nodes.dmp
```
Expected: first line `1\t|\t1\t|\tno rank\t|`.

- [ ] **Step 2: Pin the Recentrifuge container and check flags**

```bash
curl -s "https://quay.io/api/v1/repository/biocontainers/recentrifuge/tag/?limit=5&onlyActiveTags=true" | python3 -c "import sys,json;[print(t['name']) for t in json.load(sys.stdin)['tags']]"
docker run --rm quay.io/biocontainers/recentrifuge:R rcf --help | grep -E -- "-k |-c |-s |-y |-x |-e |-d "
```
Use the newest tag as `R`. Expected: all listed flags present.

- [ ] **Step 3: Write RECENTRIFUGE and its test**

`modules/local/recentrifuge/main.nf` (replace `R`):
```groovy
process RECENTRIFUGE {
    label 'process_medium'
    container 'quay.io/biocontainers/recentrifuge:R'

    input:
    val  controls     // List<String> sample ids of negative controls
    path assignments, stageAs: 'kraken/*'
    path taxdump

    output:
    path "pneuma.rcf.data.csv"   , emit: csv
    path "pneuma.rcf.html"       , emit: html
    path "pneuma_controls.rcf.*" , emit: controls_html, optional: true
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    ctrl_args=""
    other_args=""
    for f in kraken/*; do
        s=\$(basename "\$f"); s=\${s%%.*}
        case " ${controls.join(' ')} " in
            *" \$s "*) ctrl_args="\$ctrl_args -k \$f" ;;
            *)         other_args="\$other_args -k \$f" ;;
        esac
    done

    rcf -n ${taxdump} \$ctrl_args \$other_args -o pneuma.rcf.html ${args} -e CSV

    if [ ${controls.size()} -gt 0 ]; then
        rcf -n ${taxdump} \$ctrl_args \$other_args -c ${controls.size()} -o pneuma_controls.rcf.html ${args} -e CSV
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        recentrifuge: \$(rcf --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """
}
```
`modules/local/recentrifuge/tests/main.nf.test`:
```groovy
nextflow_process {
    name "Test RECENTRIFUGE"
    script "../main.nf"
    process "RECENTRIFUGE"

    test("plain + control mode") {
        setup {
            run("UNTAR") {
                script "../../../nf-core/untar/main.nf"
                process {
                    """
                    input[0] = [ [id:'db'], file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/genome/db/kraken2_bracken.tar.gz', checkIfExists: true) ]
                    """
                }
            }
            run("KRAKEN2_KRAKEN2") {
                script "../../../nf-core/kraken2/kraken2/main.nf"
                process {
                    """
                    def fq = file('https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/nanopore/fastq/test.fastq.gz', checkIfExists: true)
                    input[0] = Channel.of([ [id:'EBF0-0001'], fq ], [ [id:'CNF0-0001'], fq ])
                    input[1] = UNTAR.out.untar.map { it[1] }
                    input[2] = false
                    input[3] = true
                    """
                }
            }
        }
        when {
            process {
                """
                input[0] = ['CNF0-0001']
                input[1] = KRAKEN2_KRAKEN2.out.classified_reads_assignment.map { it[1] }.collect()
                input[2] = file("${projectDir}/tests/data/taxdump", checkIfExists: true)
                """
            }
        }
        then {
            assert process.success
            def header = path(process.out.csv[0]).readLines()[0]
            assert header.contains('EBF0-0001')
            assert process.out.controls_html.size() > 0
        }
    }
}
```
Run:
`nf-test test modules/local/recentrifuge/tests/main.nf.test --profile docker`
Expected: PASS.

- [ ] **Step 4: Write the CONTAMINATION subworkflow, params, config**

`subworkflows/local/contamination/main.nf`:
```groovy
include { RECENTRIFUGE } from '../../../modules/local/recentrifuge/main'

workflow CONTAMINATION {
    take:
    ch_assignments   // [ meta, classifiedreads.txt ]
    taxdump          // path: directory with nodes.dmp + names.dmp

    main:
    ch_controls = ch_assignments.filter { meta, f -> meta.is_control }.map { it[0].id }.collect().ifEmpty([])
    ch_files    = ch_assignments.map { it[1] }.collect()

    RECENTRIFUGE(ch_controls, ch_files, Channel.value(file(taxdump, checkIfExists: true)))

    emit:
    rcf_csv  = RECENTRIFUGE.out.csv
    html     = RECENTRIFUGE.out.html.mix(RECENTRIFUGE.out.controls_html)
    versions = RECENTRIFUGE.out.versions
}
```
`nextflow.config` params: `taxdump = null`
`conf/modules.config`:
```groovy
    withName: 'RECENTRIFUGE' {
        ext.args   = '-s NORMA -y 25 -x 9606 -d'
        publishDir = [ path: { "${params.outdir}/contamination" }, mode: params.publish_dir_mode ]
    }
```
`conf/test.config` params: `taxdump = "${projectDir}/tests/data/taxdump"`

- [ ] **Step 5: Commit**

```bash
git add tests/make_test_taxdump.py tests/data/taxdump modules/local/recentrifuge subworkflows/local/contamination nextflow.config conf
git commit -m "feat: CONTAMINATION subworkflow (Recentrifuge plain + control mode)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: FUNCTIONAL_READS subworkflow

**Files:**
- Create: `modules/local/tag_reads/main.nf`, `modules/local/abricate_reads/main.nf`, `modules/local/abricate_reads/tests/main.nf.test`, `modules/local/abricate_merge/main.nf`, `subworkflows/local/functional_reads/main.nf`
- Modify: `nextflow.config`, `conf/modules.config`

**Interfaces:**
- Consumes: `TAXONOMY.out.classified` `[ meta, classified.fastq.gz ]`.
- Produces: `FUNCTIONAL_READS.out.summary` → `path summary_reportClean.csv`; `FUNCTIONAL_READS.out.per_sample` → `[ meta, <id>_abricate.csv ]` (Phase 2 reuses it for read-level fallback); `versions`.

- [ ] **Step 1: Write the three modules**

`modules/local/tag_reads/main.nf`:
```groovy
process TAG_READS {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    tuple val(meta), path(classified)

    output:
    tuple val(meta), path("${meta.id}_kraken2ID_PFP.fasta"), emit: fasta
    path "versions.yml"                                     , emit: versions

    script:
    """
    tag_reads_taxid.py ${classified} ${meta.id}_kraken2ID_PFP.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
```
`modules/local/abricate_reads/main.nf`:
```groovy
process ABRICATE_READS {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_abricate.csv"), emit: csv
    path "versions.yml"                              , emit: versions

    script:
    def dbs = params.abricate_dbs.tokenize(',')
    """
    : > ${meta.id}_abricate.csv
    if [ -s ${fasta} ]; then
        for db in ${dbs.join(' ')}; do
            abricate ${fasta} --db \$db \\
                --minid ${params.abricate_minid} --mincov ${params.abricate_mincov} \\
                --threads ${task.cpus} --csv >> ${meta.id}_abricate.csv
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$(abricate --version | sed 's/abricate //')
        abricate_dbs: "\$(abricate --list | awk 'NR>1 {printf "%s:%s ", \$1, \$4}')"
    END_VERSIONS
    """
}
```
Verify the tag exists: `docker manifest inspect quay.io/biocontainers/abricate:1.0.1--ha8f3691_1 > /dev/null && echo ok`. If not `ok`, take the newest 1.0.1 tag from quay and use it.

`modules/local/abricate_merge/main.nf`:
```groovy
process ABRICATE_MERGE {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path csvs, stageAs: 'abricate/*'

    output:
    path "summary_reportClean.csv", emit: summary
    path "versions.yml"           , emit: versions

    script:
    """
    merge_abricate.py --out summary_reportClean.csv abricate/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Build fixtures and write the ABRICATE_READS test (empty input + real hit)**

Build a one-read FASTA that carries a CARD gene, taken from ABRicate's bundled test assembly:
```bash
IMG=quay.io/biocontainers/abricate:1.0.1--ha8f3691_1
docker run --rm -v "$PWD/tests/data:/d" $IMG bash -c '
  A=$(dirname $(readlink -f $(which abricate)))/../share/abricate*/test/assembly.fa
  [ -f $A ] || A=$(find / -name assembly.fa -path "*abricate*" 2>/dev/null | head -1)
  CONTIG=$(abricate --db card --csv $A | awk -F, "NR==2{print \$2}")
  awk -v c=">$CONTIG" "/^>/{p=(\$1==c)} p" $A \
    | sed "1s/.*/>read1_kraken:taxid|562/" > /d/S1_kraken2ID_PFP.fasta'
head -c 80 tests/data/S1_kraken2ID_PFP.fasta
: > tests/data/S0_kraken2ID_PFP.fasta
```
Expected: the first line is `>read1_kraken:taxid|562`. If no `assembly.fa` is found in the image, download `https://raw.githubusercontent.com/tseemann/abricate/master/test/assembly.fa` to `/tmp` and run the same `abricate`/`awk` lines against it.

`modules/local/abricate_reads/tests/main.nf.test`:
```groovy
nextflow_process {
    name "Test ABRICATE_READS"
    script "../main.nf"
    process "ABRICATE_READS"

    test("empty fasta gives empty csv") {
        when {
            params { abricate_dbs = 'card,vfdb,plasmidfinder'; abricate_minid = 80; abricate_mincov = 80 }
            process {
                """
                input[0] = [ [id:'S0'], file("${projectDir}/tests/data/S0_kraken2ID_PFP.fasta", checkIfExists: true) ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.csv[0][1]).text == ''
        }
    }

    test("CARD hit keeps sample file name and taxid tag") {
        when {
            params { abricate_dbs = 'card'; abricate_minid = 80; abricate_mincov = 80 }
            process {
                """
                input[0] = [ [id:'S1'], file("${projectDir}/tests/data/S1_kraken2ID_PFP.fasta", checkIfExists: true) ]
                """
            }
        }
        then {
            assert process.success
            def lines = path(process.out.csv[0][1]).readLines()
            assert lines[0].startsWith('#FILE')
            assert lines.size() > 1
            assert lines.drop(1).every { it.startsWith('S1_kraken2ID_PFP.fasta,read1_kraken:taxid|562,') }
        }
    }
}
```
Run:
`nf-test test modules/local/abricate_reads/tests/main.nf.test --profile docker`
Expected: 2 PASS.

- [ ] **Step 3: Write the subworkflow, params, config**

`subworkflows/local/functional_reads/main.nf`:
```groovy
include { TAG_READS      } from '../../../modules/local/tag_reads/main'
include { ABRICATE_READS } from '../../../modules/local/abricate_reads/main'
include { ABRICATE_MERGE } from '../../../modules/local/abricate_merge/main'

workflow FUNCTIONAL_READS {
    take:
    ch_classified   // [ meta, classified.fastq.gz ]

    main:
    TAG_READS(ch_classified)
    ABRICATE_READS(TAG_READS.out.fasta)
    ABRICATE_MERGE(ABRICATE_READS.out.csv.map { it[1] }.collect())

    emit:
    summary    = ABRICATE_MERGE.out.summary
    per_sample = ABRICATE_READS.out.csv
    versions   = TAG_READS.out.versions.first().mix(ABRICATE_READS.out.versions.first(), ABRICATE_MERGE.out.versions)
}
```
`nextflow.config` params:
```groovy
    // FUNCTIONAL_READS
    abricate_dbs    = 'card,vfdb,plasmidfinder'
    abricate_minid  = 80
    abricate_mincov = 80
```
`conf/modules.config`:
```groovy
    withName: 'ABRICATE_READS|ABRICATE_MERGE' {
        publishDir = [ path: { "${params.outdir}/functional" }, mode: params.publish_dir_mode ]
    }
    withName: 'TAG_READS' {
        publishDir = [ enabled: false ]
    }
```

- [ ] **Step 4: Commit**

```bash
git add modules/local/{tag_reads,abricate_reads,abricate_merge} subworkflows/local/functional_reads nextflow.config conf/modules.config tests/data/S0_kraken2ID_PFP.fasta tests/data/S1_kraken2ID_PFP.fasta
git commit -m "feat: FUNCTIONAL_READS subworkflow (taxid-tagged ABRicate on reads)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 10: MERIDIAN subworkflow and full pipeline wiring

**Files:**
- Create: `modules/local/meridian_config/main.nf`, `modules/local/meridian_run/main.nf`, `modules/local/meridian_run/tests/main.nf.test`, `subworkflows/local/meridian/main.nf`, `tests/main.nf.test`
- Modify: `workflows/pneuma.nf` (replace body), `main.nf` (pass DB params), `nextflow.config`, `conf/modules.config`, `conf/test.config`, `nextflow_schema.json`

**Interfaces:**
- Consumes: `TAXONOMY.out.kraken_combined`, `TAXONOMY.out.bracken_arranged`, `CONTAMINATION.out.rcf_csv`, `FUNCTIONAL_READS.out.summary`, control ids from samplesheet meta.
- Produces: `results/meridian/{Datasets,Figures,Reports}`.

- [ ] **Step 1: Write the two MERIDIAN modules**

`modules/local/meridian_config/main.nf`:
```groovy
process MERIDIAN_CONFIG {
    label 'process_single'
    container 'quay.io/biocontainers/python:3.9--1'

    input:
    path template
    path metadata
    val  controls   // List<String>

    output:
    path "config.yaml" , emit: config
    path "versions.yml", emit: versions

    script:
    def re = params.meridian_random_effect ? "--random-effect ${params.meridian_random_effect}" : ''
    """
    make_meridian_config.py \\
        --template ${template} --out config.yaml \\
        --study-id ${params.meridian_study_id} \\
        --metadata ${metadata} \\
        --group ${params.meridian_group} ${re} \\
        --controls "${controls.join(',')}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
```
`modules/local/meridian_run/main.nf`:
```groovy
process MERIDIAN_RUN {
    label 'process_high'
    container 'ghcr.io/julio92-c/meridian:v1.0.0'
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
    ( cd /opt/meridian && Rscript run_pipeline.R "\$OLDPWD/run_config.yaml" )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meridian: v1.0.0
    END_VERSIONS
    """

    stub:
    """
    mkdir -p Datasets Figures Reports
    touch Reports/pipeline_report.html
    echo '"${task.process}": {meridian: v1.0.0}' > versions.yml
    """
}
```
`modules/local/meridian_run/tests/main.nf.test`:
```groovy
nextflow_process {
    name "Test MERIDIAN_RUN (stub)"
    script "../main.nf"
    process "MERIDIAN_RUN"
    options "-stub"

    test("stages inputs and produces report dirs") {
        when {
            process {
                """
                input[0] = file("${projectDir}/assets/meridian_config.template.yaml")
                input[1] = file("${projectDir}/tests/data/metadata_test.csv")
                input[2] = file("${projectDir}/tests/data/abricate_S1.csv")
                input[3] = file("${projectDir}/tests/data/kraken_combined.txt")
                input[4] = file("${projectDir}/tests/data/bracken_combined.txt")
                input[5] = file("${projectDir}/tests/data/abricate_S2_empty.csv")
                """
            }
        }
        then {
            assert process.success
            assert path("${process.out.reports[0]}/pipeline_report.html").exists()
        }
    }
}
```
Run: `nf-test test modules/local/meridian_run/tests/main.nf.test --profile docker`
Expected: PASS.

- [ ] **Step 2: Write the MERIDIAN subworkflow**

`subworkflows/local/meridian/main.nf`:
```groovy
include { MERIDIAN_CONFIG } from '../../../modules/local/meridian_config/main'
include { MERIDIAN_RUN    } from '../../../modules/local/meridian_run/main'

workflow MERIDIAN {
    take:
    ch_controls        // value channel: List<String>
    abricate_summary
    kraken_combined
    bracken_arranged
    rcf_csv

    main:
    def metadata = file(params.meridian_metadata, checkIfExists: true)
    MERIDIAN_CONFIG(file("${projectDir}/assets/meridian_config.template.yaml"), metadata, ch_controls)
    MERIDIAN_RUN(MERIDIAN_CONFIG.out.config, metadata, abricate_summary, kraken_combined, bracken_arranged, rcf_csv)

    emit:
    reports  = MERIDIAN_RUN.out.reports
    versions = MERIDIAN_CONFIG.out.versions.mix(MERIDIAN_RUN.out.versions)
}
```

- [ ] **Step 3: Wire the pipeline**

Replace the body of `workflows/pneuma.nf` (keep the template's MultiQC / versions blocks at the end of `main:`; replace everything above them) with:
```groovy
include { PREPROCESS       } from '../subworkflows/local/preprocess/main'
include { TAXONOMY         } from '../subworkflows/local/taxonomy/main'
include { CONTAMINATION    } from '../subworkflows/local/contamination/main'
include { FUNCTIONAL_READS } from '../subworkflows/local/functional_reads/main'
include { MERIDIAN         } from '../subworkflows/local/meridian/main'

workflow PNEUMA {
    take:
    ch_samplesheet   // [ meta, [fastq...] ]

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    PREPROCESS(ch_samplesheet, file(params.host_index, checkIfExists: true))
    TAXONOMY(PREPROCESS.out.reads, params.kraken2_db)
    CONTAMINATION(TAXONOMY.out.assignments, params.taxdump)
    FUNCTIONAL_READS(TAXONOMY.out.classified)

    ch_versions = ch_versions.mix(PREPROCESS.out.versions, TAXONOMY.out.versions,
                                  CONTAMINATION.out.versions, FUNCTIONAL_READS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS.out.multiqc)

    if (!params.skip_meridian) {
        ch_controls = ch_samplesheet.filter { it[0].is_control }.map { it[0].id }.collect().ifEmpty([])
        MERIDIAN(ch_controls, FUNCTIONAL_READS.out.summary, TAXONOMY.out.kraken_combined,
                 TAXONOMY.out.bracken_arranged, CONTAMINATION.out.rcf_csv)
        ch_versions = ch_versions.mix(MERIDIAN.out.versions)
    }

    if (params.database_manifest) {
        Channel.fromPath(params.database_manifest, checkIfExists: true)
            .collectFile(name: 'database_manifest.tsv', storeDir: "${params.outdir}/pipeline_info")
    }

    // ---- template: softwareVersionsToYAML + MULTIQC block continues here ----
```
Add to `nextflow.config` params:
```groovy
    // MERIDIAN
    skip_meridian          = false
    meridian_metadata      = null
    meridian_group         = null
    meridian_random_effect = null
    meridian_study_id      = 'pneuma'
    // Provenance
    database_manifest      = null
```
Add to `conf/modules.config`:
```groovy
    withName: 'MERIDIAN_RUN' {
        publishDir = [ path: { "${params.outdir}/meridian" }, mode: params.publish_dir_mode ]
    }
    withName: 'MERIDIAN_CONFIG' {
        publishDir = [ path: { "${params.outdir}/meridian" }, mode: params.publish_dir_mode ]
    }
```
Add to `conf/test.config` params: `skip_meridian = true` (3 test samples are too few for MERIDIAN's statistics; MERIDIAN integration is validated on real data in Task 12).

Add a parameter guard to `PIPELINE_INITIALISATION` (after `validateParameters()`):
```groovy
    ['host_index', 'kraken2_db', 'taxdump'].each { p ->
        if (!params[p]) error("--${p} is required")
    }
    if (!params.skip_meridian && !(params.meridian_metadata && params.meridian_group)) {
        error("--meridian_metadata and --meridian_group are required unless --skip_meridian is set")
    }
```

- [ ] **Step 4: Update the parameter schema**

Run: `nf-core pipelines schema build`
Answer `y` to add each new parameter, `n` to launching the web builder. Then edit `nextflow_schema.json` so `host_index`, `kraken2_db`, `taxdump`, `meridian_metadata` have `"format": "path"` (or `"file-path"`), and every new param has a one-line `description`.
Run: `nf-core pipelines lint --dir . -k nextflow_config -k schema_params`
Expected: no failures for those checks.

- [ ] **Step 5: Write and run the pipeline-level test**

`tests/main.nf.test`:
```groovy
nextflow_pipeline {
    name "PNEUMA end-to-end (test profile)"
    script "../main.nf"
    profile "test"

    test("core outputs exist") {
        when { params { outdir = "$outputDir" } }
        then {
            assert workflow.success
            assert path("$outputDir/taxonomy/bracken_arranged.csv").readLines()[0] ==
                   '#perc,tot_all,CNF0-0001,EBF0-0001,SPF0-0001,taxid,name'
            assert path("$outputDir/taxonomy/kraken2_combined_reports.txt").exists()
            assert path("$outputDir/contamination/pneuma.rcf.data.csv").exists()
            assert path("$outputDir/functional/summary_reportClean.csv").readLines()[0].startsWith('#FILE')
            assert path("$outputDir/multiqc/multiqc_report.html").exists()
        }
    }
}
```
Run: `nf-test test tests/main.nf.test --profile docker`
Expected: PASS. Then run the whole suite: `nf-test test --profile docker && pytest tests/python -v` — all PASS.

- [ ] **Step 6: Commit**

```bash
git add -A
git commit -m "feat: MERIDIAN handoff and end-to-end pipeline wiring

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 11: Myriad profile, database preparation, docs

**Files:**
- Create: `conf/myriad.config`, `bin/prepare_databases.sh`, `tests/python/test_prepare_databases.py`, `docs/usage.md` (replace template text)
- Modify: `nextflow.config` (profile `myriad`)

**Interfaces:**
- Produces: `bin/prepare_databases.sh --outdir DIR --kraken-release YYYYMMDD --taxdump-month YYYY-MM [--dry-run]` writing `DIR/database_manifest.tsv` (columns `name\tversion\tpath\tsha256`) and printing the four `--param` lines to pass to the pipeline.

- [ ] **Step 1: Write the failing test for the dry run**

`tests/python/test_prepare_databases.py`:
```python
import subprocess
from conftest import BIN


def run(*args):
    return subprocess.run(["bash", str(BIN / "prepare_databases.sh"), *args],
                          capture_output=True, text=True)


def test_dry_run_prints_pinned_urls(tmp_path):
    r = run("--outdir", str(tmp_path), "--kraken-release", "20250402",
            "--taxdump-month", "2025-04", "--dry-run")
    assert r.returncode == 0, r.stderr
    assert "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20250402.tar.gz" in r.stdout
    assert "taxdmp_2025-04-01.zip" in r.stdout
    assert "hostile fetch --aligner minimap2 --name human-t2t-hla" in r.stdout


def test_missing_args_fail(tmp_path):
    r = run("--outdir", str(tmp_path))
    assert r.returncode != 0
    assert "--kraken-release" in r.stderr
```
Run: `pytest tests/python/test_prepare_databases.py -v` → FAIL (script missing).

- [ ] **Step 2: Implement `bin/prepare_databases.sh`**

```bash
#!/usr/bin/env bash
# Download and pin PNEUMA Phase-1 databases; write database_manifest.tsv.
# Usage: prepare_databases.sh --outdir DIR --kraken-release YYYYMMDD --taxdump-month YYYY-MM
#                             [--hostile-image IMG] [--dry-run]
set -euo pipefail

OUTDIR="" ; KRELEASE="" ; TMONTH="" ; DRY=0
HOSTILE_IMAGE="docker://quay.io/biocontainers/hostile:T"   # same tag as modules/local/hostile_clean

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)         OUTDIR="$2"; shift 2 ;;
    --kraken-release) KRELEASE="$2"; shift 2 ;;
    --taxdump-month)  TMONTH="$2"; shift 2 ;;
    --hostile-image)  HOSTILE_IMAGE="$2"; shift 2 ;;
    --dry-run)        DRY=1; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$OUTDIR"   ]] || { echo "--outdir is required" >&2; exit 2; }
[[ -n "$KRELEASE" ]] || { echo "--kraken-release YYYYMMDD is required (see https://benlangmead.github.io/aws-indexes/k2)" >&2; exit 2; }
[[ -n "$TMONTH"   ]] || { echo "--taxdump-month YYYY-MM is required (match the Kraken2 release month)" >&2; exit 2; }

K_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_${KRELEASE}.tar.gz"
T_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${TMONTH}-01.zip"

run() { if [[ $DRY -eq 1 ]]; then echo "$*"; else eval "$@"; fi; }
sha() { if [[ $DRY -eq 1 ]]; then echo "dry-run"; else sha256sum "$1" | cut -d' ' -f1; fi; }

run "mkdir -p '$OUTDIR/kraken2_pluspf_${KRELEASE}' '$OUTDIR/taxdump_${TMONTH}' '$OUTDIR/hostile'"
run "wget -nv -O '$OUTDIR/k2_pluspf_${KRELEASE}.tar.gz' '$K_URL'"
run "tar -xzf '$OUTDIR/k2_pluspf_${KRELEASE}.tar.gz' -C '$OUTDIR/kraken2_pluspf_${KRELEASE}'"
run "wget -nv -O '$OUTDIR/taxdmp_${TMONTH}-01.zip' '$T_URL'"
run "unzip -o -q '$OUTDIR/taxdmp_${TMONTH}-01.zip' nodes.dmp names.dmp -d '$OUTDIR/taxdump_${TMONTH}'"
run "HOSTILE_CACHE_DIR='$OUTDIR/hostile' singularity exec $HOSTILE_IMAGE hostile fetch --aligner minimap2 --name human-t2t-hla"

MANIFEST="$OUTDIR/database_manifest.tsv"
if [[ $DRY -eq 0 ]]; then
  HOSTILE_MMI=$(find "$OUTDIR/hostile" -name 'human-t2t-hla*.mmi' | head -1)
  {
    printf 'name\tversion\tpath\tsha256\n'
    printf 'kraken2_pluspf\t%s\t%s\t%s\n' "$KRELEASE" "$OUTDIR/kraken2_pluspf_${KRELEASE}" "$(sha "$OUTDIR/k2_pluspf_${KRELEASE}.tar.gz")"
    printf 'ncbi_taxdump\t%s\t%s\t%s\n' "$TMONTH" "$OUTDIR/taxdump_${TMONTH}" "$(sha "$OUTDIR/taxdmp_${TMONTH}-01.zip")"
    printf 'hostile_index\thuman-t2t-hla\t%s\t%s\n' "$HOSTILE_MMI" "$(sha "$HOSTILE_MMI")"
  } > "$MANIFEST"
  echo "Wrote $MANIFEST"
  echo "--kraken2_db $OUTDIR/kraken2_pluspf_${KRELEASE}"
  echo "--taxdump $OUTDIR/taxdump_${TMONTH}"
  echo "--host_index $HOSTILE_MMI"
  echo "--database_manifest $MANIFEST"
fi
```
Replace `T` with the Hostile tag from Task 6. `chmod +x bin/prepare_databases.sh`.
Run: `pytest tests/python -v` → all PASS.

- [ ] **Step 3: Write `conf/myriad.config` and register the profile**

`conf/myriad.config`:
```groovy
params {
    config_profile_name        = 'UCL Myriad (SGE)'
    config_profile_description = 'PNEUMA on UCL Myriad: SGE, Singularity, per-core memory requests'
}

executor {
    name            = 'sge'
    queueSize       = 50
    submitRateLimit = '10 sec'
}

process {
    executor = 'sge'
    penv     = 'smp'
    // Myriad's -l mem is PER CORE. Nextflow adds -pe smp <cpus> and -l h_rt from task.time.
    clusterOptions = { "-l mem=${Math.max(1, Math.ceil(task.memory.toGiga() / task.cpus) as int)}G -l tmpfs=20G" }
    resourceLimits = [ cpus: 36, memory: 1500.GB, time: 48.h ]

    withName: 'KRAKEN2_KRAKEN2' {
        cpus   = 12
        memory = 110.GB
        time   = 24.h
    }
}

singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = "${System.getenv('HOME')}/Scratch/singularity_cache"
}
```
In `nextflow.config` `profiles { }` add: `myriad { includeConfig 'conf/myriad.config' }`

- [ ] **Step 4: Verify the SGE headers Nextflow generates (on Myriad login node)**

```bash
nextflow run Julio92-C/pneuma -r main -profile myriad,test --outdir /tmp/pneuma_stub -stub-run
grep '^#\$' work/*/*/.command.run | sort -u | head -20
```
Expected: lines containing `-pe smp`, `-l h_rt=`, `-l mem=<n>G`, `-l tmpfs=20G`.
Then submit one real job: `nextflow run ~/Scratch/pneuma -profile myriad,test --outdir /tmp/pneuma_t -resume`.
If `qsub` rejects a memory resource that Nextflow added itself (e.g. `h_rss` or `virtual_free`), record the exact rejected header line in `docs/usage.md` under "Known Myriad issues" and stop: that is a design decision for Julio (options: per-process `clusterOptions` only, or asking Myriad support), not something to work around silently.

- [ ] **Step 5: Write `docs/usage.md`**

```markdown
# PNEUMA usage

## 1. Stage databases (once, on Myriad login node)
    bin/prepare_databases.sh --outdir ~/Scratch/pneuma_db --kraken-release <YYYYMMDD> --taxdump-month <YYYY-MM>
Copy the four `--param` lines it prints.

## 2. Samplesheet
`sample,run_id,batch,site,fastq,is_control` — one row per sample × run. A sample sequenced on
several flow cells gets several rows with the same `sample`; their reads are concatenated.
`sample` must match `^[A-Za-z0-9][A-Za-z0-9-]*$`.

## 3. Run
    nextflow run Julio92-C/pneuma -r <tag> -profile myriad \
      --input samplesheet.csv --outdir results \
      --kraken2_db ... --taxdump ... --host_index ... --database_manifest ... \
      --meridian_metadata metadata.csv --meridian_group sample_type --meridian_random_effect ID \
      -resume

## 4. Outputs
- `preprocess/` nanoq stats, host-free reads, NanoPlot, seqkit stats
- `taxonomy/` Kraken2 reports + per-read assignments, Bracken, `bracken_arranged.csv`, `kraken2_combined_reports.txt`, Krona
- `contamination/` `pneuma.rcf.data.csv`, `pneuma.rcf.html`, `pneuma_controls.rcf.html`
- `functional/` per-sample ABRicate, `summary_reportClean.csv`
- `meridian/` MERIDIAN Datasets/Figures/Reports
- `pipeline_info/` software versions, `database_manifest.tsv`
```

- [ ] **Step 6: Commit**

```bash
git add conf/myriad.config bin/prepare_databases.sh tests/python/test_prepare_databases.py docs/usage.md nextflow.config
git commit -m "feat: Myriad profile, database preparation script, usage docs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 12: Myriad smoke run on real data (validation gate)

**Files:**
- Create: `docs/validation/phase1_smoke.md` (results write-up). No code changes unless a defect is found (fix via its own TDD cycle, then re-run).

- [ ] **Step 1: Stage databases on Myriad**

```bash
ssh myriad
cd ~/Scratch && git clone git@github.com:Julio92-C/pneuma.git && cd pneuma
module load java/temurin-17 singularity-env   # or the module names Myriad lists via `module avail`
bin/prepare_databases.sh --outdir ~/Scratch/pneuma_db --kraken-release <latest PlusPF date> --taxdump-month <same YYYY-MM>
```
Expected: `database_manifest.tsv` with three rows and real sha256 values.

- [ ] **Step 2: Build a 6-sample smoke samplesheet from EBSP 2026**

Use two RFH runs: 2 EB + 1 SP from one run, 1 EB + 1 control (`CNN0-0609`) from the last run, and one sample sequenced on two flow cells (two rows, same `sample`), pointing at `~/Scratch/JC_2026-09_EBSP/merged_fastq/<sample>/*.fastq.gz`. Save as `~/Scratch/pneuma_smoke/samplesheet.csv`.

- [ ] **Step 3: Run**

```bash
cd ~/Scratch/pneuma_smoke
nextflow run ~/Scratch/pneuma -profile myriad --input samplesheet.csv --outdir results \
  <four --param lines from prepare_databases.sh> --skip_meridian -resume
```
Expected: `Pipeline completed successfully`; warning lists runs without controls.

- [ ] **Step 4: Check against legacy outputs for the same samples**

For each smoke sample, compare `results/taxonomy/bracken_arranged.csv` with `JC_2026-09_EBSP/Reports/pipeline_inputs/bracken_arranged.csv`: top-10 species by count overlap, and total classified reads. Compare `results/functional/summary_reportClean.csv` gene sets with the legacy `summary_reportClean.csv` per sample. Record the numbers in `docs/validation/phase1_smoke.md`, noting that differences are expected from PlusPF-8 → full PlusPF and from Hostile vs. minimap2/GRCh38.

- [ ] **Step 5: MERIDIAN on the full RFH 2026 cohort**

Build the full EBSP 2026 samplesheet (89 analysis samples + control, re-extraction EBR0-2466 given `sample` = `EBF0-2466`) and run with MERIDIAN on:
```bash
nextflow run ~/Scratch/pneuma -profile myriad --input samplesheet_rfh2026.csv --outdir results_rfh2026 \
  <four --param lines> --meridian_metadata sample_metadata.csv --meridian_group sample_type \
  --meridian_random_effect ID -resume
```
Expected: `results_rfh2026/meridian/Reports/pipeline_report.html` renders. Note in `phase1_smoke.md` any MERIDIAN stage that errors, with the log excerpt.

- [ ] **Step 6: Commit the write-up and tag**

```bash
git add docs/validation/phase1_smoke.md
git commit -m "docs: Phase 1 Myriad smoke validation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git tag v0.1.0 && git push --follow-tags
```
