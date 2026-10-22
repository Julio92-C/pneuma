#!/usr/bin/env python3
"""Write one database_manifest row (name, version, path, sha256) for a prepared database.

  database_row.py --name bakta_db --rel bakta/db-light --target db-light \
      (--version 6.0 | --version-from bakta|amrfinderplus) --out 08_bakta_db.row.tsv

sha256 is the file digest, or for a directory the digest of the sorted
'<sha256>  <relative path>' lines of all its files (stable across copies and timestamps).
"""
import argparse
import hashlib
import json
import os
import sys
import tarfile


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sha256_path(path):
    path = str(path)
    if os.path.isfile(path):
        return sha256_file(path)
    lines = []
    for root, _dirs, files in os.walk(path):
        for name in files:
            full = os.path.join(root, name)
            lines.append("%s  %s" % (sha256_file(full), os.path.relpath(full, path).replace(os.sep, "/")))
    return hashlib.sha256(("\n".join(sorted(lines)) + "\n").encode()).hexdigest()


def version_from(kind, target):
    if kind == "bakta":
        with open(os.path.join(str(target), "version.json")) as fh:
            v = json.load(fh)
        return "%s.%s-%s (%s)" % (v["major"], v["minor"], v["type"], v["date"])
    if kind == "amrfinderplus":
        with tarfile.open(str(target), "r:gz") as t:
            member = next(m for m in t.getmembers() if os.path.basename(m.name.rstrip("/")) == "version.txt")
            return t.extractfile(member).read().decode().strip()
    raise ValueError("unknown --version-from %s" % kind)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--rel", required=True, help="path relative to --db_outdir, as written in the manifest")
    ap.add_argument("--target", required=True, help="file or directory to checksum")
    ap.add_argument("--out", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--version")
    g.add_argument("--version-from", choices=["bakta", "amrfinderplus"])
    args = ap.parse_args(argv)
    version = args.version if args.version is not None else version_from(args.version_from, args.target)
    if not version or not version.strip():
        ap.error("empty version for %s" % args.name)
    sha = sha256_path(args.target)
    with open(args.out, "w") as out:
        out.write("%s\t%s\t%s\t%s\n" % (args.name, version.strip(), args.rel, sha))
    return 0


if __name__ == "__main__":
    sys.exit(main())
