#!/usr/bin/env python3
"""Render a MERIDIAN v1.0.0 config.yaml from PNEUMA's template.

Placeholders: {{STUDY_ID}} {{STUDY_NAME}} {{METADATA_FILE}} {{GROUP_COL}} {{RANDOM_EFFECT}} {{CONTROLS}}.
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


def render(template_text, study_id, metadata_name, group, random_effect, controls):
    text = (template_text
            .replace("{{STUDY_ID}}", yaml_str(study_id))
            .replace("{{STUDY_NAME}}", yaml_str(f"{study_id} (PNEUMA)"))
            .replace("{{METADATA_FILE}}", yaml_str(f"inputs/{metadata_name}"))
            .replace("{{GROUP_COL}}", yaml_str(group))
            .replace("{{RANDOM_EFFECT}}", yaml_str(random_effect) if random_effect else "~")
            .replace("{{CONTROLS}}", yaml_list(controls)))
    if "{{" in text:
        raise ValueError("unfilled placeholder left in MERIDIAN config template")
    return text


def check_metadata(path, group, random_effect, controls):
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
    args = parser.parse_args(argv)

    controls = [c.strip() for c in args.controls.split(",") if c.strip()]
    check_metadata(args.metadata, args.group, args.random_effect, controls)
    text = render(Path(args.template).read_text(encoding="utf-8"), args.study_id,
                  Path(args.metadata).name, args.group, args.random_effect, controls)
    Path(args.out).write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
