#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample", required=True)
    return parser.parse_args()


def open_text(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def parse_info(info_field):
    info_items = []
    info_map = {}
    if info_field == ".":
        return info_items, info_map

    for item in info_field.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        else:
            key, value = item, None
        info_items.append((key, value))
        info_map[key] = value
    return info_items, info_map


def format_info(info_items):
    if not info_items:
        return "."
    parts = []
    for key, value in info_items:
        parts.append(f"{key}={value}" if value is not None else key)
    return ";".join(parts)


def add_missing_info(info_items, info_map):
    if "DP4" in info_map and info_map["DP4"]:
        try:
            dp4_values = [int(x) for x in info_map["DP4"].split(",")]
        except ValueError:
            dp4_values = []
    else:
        dp4_values = []

    dp_value = info_map.get("DP")
    if (dp_value is None or dp_value == "") and dp4_values:
        dp_value = str(sum(dp4_values))
        info_items.append(("DP", dp_value))
        info_map["DP"] = dp_value

    af_value = info_map.get("AF")
    if (af_value is None or af_value == "") and dp4_values:
        depth = sum(dp4_values)
        alt_depth = dp4_values[2] + dp4_values[3] if len(dp4_values) >= 4 else 0
        af = (alt_depth / depth) if depth else 0.0
        af_value = f"{af:.6f}".rstrip("0").rstrip(".")
        if not af_value:
            af_value = "0"
        info_items.append(("AF", af_value))
        info_map["AF"] = af_value

    return info_items


def main():
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    seen_format_gt = False
    seen_info_dp = False
    seen_info_af = False

    with open_text(input_path, "r") as in_handle, open(output_path, "w", encoding="utf-8") as out_handle:
        for raw_line in in_handle:
            line = raw_line.rstrip("\n")

            if line.startswith("##FORMAT=<ID=GT,"):
                seen_format_gt = True
                out_handle.write(raw_line)
                continue

            if line.startswith("##INFO=<ID=DP,"):
                seen_info_dp = True
                out_handle.write(raw_line)
                continue

            if line.startswith("##INFO=<ID=AF,"):
                seen_info_af = True
                out_handle.write(raw_line)
                continue

            if line.startswith("#CHROM"):
                if not seen_info_dp:
                    out_handle.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
                if not seen_info_af:
                    out_handle.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
                if not seen_format_gt:
                    out_handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

                columns = line.split("\t")
                if len(columns) == 8:
                    columns.extend(["FORMAT", args.sample])
                out_handle.write("\t".join(columns) + "\n")
                continue

            if line.startswith("#"):
                out_handle.write(raw_line)
                continue

            fields = line.split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line}")

            info_items, info_map = parse_info(fields[7])
            fields[7] = format_info(add_missing_info(info_items, info_map))

            if len(fields) == 8:
                fields.extend(["GT", "1"])
            else:
                if len(fields) < 10:
                    raise ValueError(f"Unexpected VCF column count: {line}")
                fields[8] = "GT"
                fields[9] = "1"

            out_handle.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
