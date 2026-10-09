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
