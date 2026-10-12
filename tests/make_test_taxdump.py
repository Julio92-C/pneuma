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
        nodes.append(f"{taxid}\t|\t{parent if taxid != '1' else '1'}\t|\t{rank}\t|\t\t|")
        names.append(f"{taxid}\t|\t{raw.strip()}\t|\t\t|\tscientific name\t|")
        stack.append((depth, taxid))
    # rcf's -x 9606 (exclude human) requires the taxid to exist in the dump
    if not any(n.startswith("9606	") for n in nodes):
        nodes.append("9606	|	1	|	species	|		|")
        names.append("9606	|	Homo sapiens	|		|	scientific name	|")
    (out / "nodes.dmp").write_text("\n".join(nodes) + "\n")
    (out / "names.dmp").write_text("\n".join(names) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
