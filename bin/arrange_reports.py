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


def sample_names_from_preamble(lines):
    """Sample names written by `combine_kreports.py --sample-names`: '#<name>\t<path>' lines."""
    names = []
    for line in lines:
        if line.startswith("#perc\t"):
            break
        if line.startswith("#") and "\t" in line:
            names.append(line[1:].split("\t")[0])
    return names


def fmt_number(value):
    value = round(value, 4)
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.4f}".rstrip("0").rstrip(".")


def arrange_bracken(lines, names=None):
    header = lines[0].split("\t")
    sample_cols = [c for c in header if c.endswith("_all") and c != "tot_all"]
    samples = [c[: -len("_all")] for c in sample_cols]
    # combine_kreports.py labels the columns 1_all, 2_all, ...; real names are in the preamble
    if names and len(names) == len(samples) and all(x.isdigit() for x in samples):
        samples = list(names)
    idx = {c: header.index(c) for c in sample_cols + ["taxid", "name"]}
    body = []
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < len(header) and fields[-3:-2] == ["U"]:
            # combine_kreports.py writes a totals-only unclassified row when no report has one
            fields = fields[:3] + ["0"] * (len(header) - 6) + fields[-3:]
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
        raw = [line.rstrip("\r\n") for line in fh.readlines()]
    rows = arrange_bracken(read_from_perc(raw), sample_names_from_preamble(raw))
    with open(args.bracken_out, "w", newline="") as out:
        csv.writer(out, lineterminator="\n").writerows(rows)


if __name__ == "__main__":
    main()
