#!/usr/bin/env python3
"""Render a MERIDIAN (>= 1.2.0) config.yaml from PNEUMA's template.

Placeholders: {{STUDY_ID}} {{STUDY_NAME}} {{METADATA_FILE}} {{GROUP_COL}} {{RANDOM_EFFECT}} {{CONTROLS}}
{{CONTROL_MATCH_COLS}} {{RANDOM_EFFECT_TESTS}} {{MIN_SAMPLE_READS}}.
String values are emitted as double-quoted YAML scalars so values such as
yes/no/null/2024 stay strings. '__PROJECT_ROOT__' is left for the MERIDIAN
process to replace with its work dir.
"""
import argparse
import csv
from pathlib import Path


NA_VALUES = {"", "na", "n/a", "nan", "null", "none"}


def yaml_str(value):
    escaped = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def yaml_list(items):
    return "[" + ", ".join(yaml_str(i) for i in items) + "]"


def render(template_text, study_id, metadata_name, group, random_effect, controls,
           control_match_cols=(), random_effect_tests=False, min_sample_reads=0):
    text = (template_text
            .replace("{{STUDY_ID}}", yaml_str(study_id))
            .replace("{{STUDY_NAME}}", yaml_str(f"{study_id} (PNEUMA)"))
            .replace("{{METADATA_FILE}}", yaml_str(f"inputs/{metadata_name}"))
            .replace("{{GROUP_COL}}", yaml_str(group))
            .replace("{{RANDOM_EFFECT}}", yaml_str(random_effect) if random_effect else "~")
            .replace("{{CONTROLS}}", yaml_list(controls))
            .replace("{{CONTROL_MATCH_COLS}}", yaml_list(control_match_cols))
            .replace("{{RANDOM_EFFECT_TESTS}}", "true" if random_effect_tests else "false")
            .replace("{{MIN_SAMPLE_READS}}", str(int(min_sample_reads))))
    if "{{" in text:
        raise ValueError("unfilled placeholder left in MERIDIAN config template")
    return text


def check_metadata(path, group, random_effect, controls, control_match_cols=()):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        columns = reader.fieldnames or []
        rows = list(reader) if "sample" in columns else []
        samples = {row["sample"] for row in rows}
    if "sample" not in columns:
        raise ValueError(f"metadata {path} has no 'sample' column")
    if group not in columns:
        raise ValueError(f"group column '{group}' not in metadata columns {columns}")
    if random_effect and random_effect not in columns:
        raise ValueError(f"random-effect column '{random_effect}' not in metadata columns {columns}")
    for col in control_match_cols:
        if col not in columns:
            raise ValueError(f"control-match column '{col}' not in metadata columns {columns}")
    missing = [c for c in controls if c not in samples]
    if missing:
        raise ValueError(f"control sample(s) not in metadata: {missing}")
    levels = {(row[group] or "").strip() for row in rows if row["sample"] not in controls}
    levels = {v for v in levels if v.lower() not in NA_VALUES}
    if len(levels) < 2:
        raise ValueError(
            f"group column '{group}' needs at least 2 distinct non-NA levels among the non-control "
            f"samples, found {len(levels)} ({sorted(levels)})")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--group", required=True)
    parser.add_argument("--random-effect", default=None)
    parser.add_argument("--controls", default="")
    parser.add_argument("--control-match-cols", default="",
                        help="comma-separated metadata columns; each sample uses the controls sharing its value "
                             "in the first column that has any (e.g. seq_run,batch)")
    parser.add_argument("--random-effect-tests", action="store_true",
                        help="patient/cluster-aware permutation tests on --random-effect")
    parser.add_argument("--min-sample-reads", type=int, default=0,
                        help="beta diversity: drop samples with fewer reads (0 = keep all)")
    args = parser.parse_args(argv)

    controls = [c.strip() for c in args.controls.split(",") if c.strip()]
    match_cols = [c.strip() for c in args.control_match_cols.split(",") if c.strip()]
    if args.random_effect_tests and not args.random_effect:
        raise ValueError("random_effect_tests needs --random-effect (the column to permute within)")
    if args.min_sample_reads < 0:
        raise ValueError("--min-sample-reads must be >= 0")
    check_metadata(args.metadata, args.group, args.random_effect, controls, match_cols)
    text = render(Path(args.template).read_text(encoding="utf-8"), args.study_id,
                  Path(args.metadata).name, args.group, args.random_effect, controls,
                  match_cols, args.random_effect_tests, args.min_sample_reads)
    Path(args.out).write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
