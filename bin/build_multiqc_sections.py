#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Build custom MultiQC sections from pipeline outputs.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--coverage-output", required=True)
    parser.add_argument("--read-output", required=True)
    parser.add_argument("--variant-output", required=True)
    return parser.parse_args()


def infer_sample_and_read_type(path):
    parts = Path(path).parts
    sample_id = None
    read_type = None
    for i, part in enumerate(parts):
        if part in {"long_reads", "short_reads"} and i > 0:
            sample_id = parts[i - 1]
            read_type = "long" if part == "long_reads" else "short"
            break
    if sample_id is None:
        name = Path(path).name
        sample_id = name.split(".")[0]
        read_type = "unknown"
    return sample_id, read_type


def row_label(sample_id, read_type):
    if read_type in {"long", "short"}:
        return f"{sample_id} ({read_type})"
    return sample_id


def load_tsv_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_seqkit_row(path):
    rows = load_tsv_rows(path)
    if not rows:
        raise ValueError(f"No seqkit rows found in {path}")
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one seqkit row in {path}, found {len(rows)}")
    return rows[0]


def to_float(value, default=None):
    if value in {"", None}:
        return default
    return float(value)


def to_int(value, default=None):
    if value in {"", None}:
        return default
    return int(float(value))


def parse_mosdepth_global_dist(path):
    target_fraction = None
    median_depth = 0

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            chrom, depth, fraction = line.split("\t")
            if chrom != "total":
                continue
            depth = int(depth)
            fraction = float(fraction)
            if depth == 10:
                target_fraction = fraction * 100.0
            if fraction >= 0.5 and depth > median_depth:
                median_depth = depth

    return {
        "coverage_10x_pct": target_fraction if target_fraction is not None else 0.0,
        "median_coverage": median_depth,
    }


def parse_mosdepth_summary(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["chrom"] == "total":
                return {
                    "mean_coverage": to_float(row.get("mean"), 0.0),
                    "min_coverage": to_float(row.get("min"), 0.0),
                    "max_coverage": to_float(row.get("max"), 0.0),
                    "genome_length": to_int(row.get("length"), 0),
                }
    raise ValueError(f"No total row found in {path}")


def build_coverage_section(outdir):
    rows = {}

    for path in sorted(outdir.rglob("*.mosdepth.summary.txt")):
        sample_id, read_type = infer_sample_and_read_type(path)
        key = (sample_id, read_type)
        rows.setdefault(key, {"sample_id": sample_id, "read_type": read_type})
        rows[key].update(parse_mosdepth_summary(path))

    for path in sorted(outdir.rglob("*.mosdepth.global.dist.txt")):
        sample_id, read_type = infer_sample_and_read_type(path)
        key = (sample_id, read_type)
        rows.setdefault(key, {"sample_id": sample_id, "read_type": read_type})
        rows[key].update(parse_mosdepth_global_dist(path))

    table_data = {}
    for key in sorted(rows):
        row = rows[key]
        table_data[row_label(row["sample_id"], row["read_type"])] = {
            "coverage_10x_pct": row.get("coverage_10x_pct", ""),
            "median_coverage": row.get("median_coverage", ""),
            "mean_coverage": row.get("mean_coverage", ""),
            "min_coverage": row.get("min_coverage", ""),
            "max_coverage": row.get("max_coverage", ""),
            "genome_length": row.get("genome_length", ""),
        }

    return {
        "id": "h2seq_coverage_statistics",
        "section_name": "Coverage Statistics",
        "description": "Coverage summary metrics derived from mosdepth outputs.",
        "plot_type": "table",
        "pconfig": {
            "id": "h2seq_coverage_statistics_table",
            "title": "Coverage Statistics",
        },
        "headers": {
            "coverage_10x_pct": {"title": "≥ 10X", "format": "{:,.1f}", "suffix": "%"},
            "median_coverage": {"title": "Median", "format": "{:,.1f}X"},
            "mean_coverage": {"title": "Mean Cov.", "format": "{:,.1f}"},
            "min_coverage": {"title": "Min Cov.", "format": "{:,.1f}"},
            "max_coverage": {"title": "Max Cov.", "format": "{:,.1f}"},
            "genome_length": {"title": "Genome length"},
        },
        "data": table_data,
    }


def build_read_section(outdir):
    rows = {}

    for path in sorted(outdir.rglob("*.raw_long.tsv")) + sorted(outdir.rglob("*.raw_short.tsv")):
        sample_id, read_type = infer_sample_and_read_type(path)
        record = load_seqkit_row(path)
        rows[(sample_id, read_type)] = {
            "sample_id": sample_id,
            "read_type": read_type,
            "raw_reads": to_int(record.get("num_seqs"), ""),
            "raw_bases": to_int(record.get("sum_len"), ""),
            "raw_mean_length": to_float(record.get("avg_len"), ""),
            "raw_mean_quality": to_float(record.get("AvgQual"), ""),
        }

    for path in sorted(outdir.rglob("*.clean_long.tsv")) + sorted(outdir.rglob("*.clean_short.tsv")):
        sample_id, read_type = infer_sample_and_read_type(path)
        record = load_seqkit_row(path)
        rows.setdefault((sample_id, read_type), {"sample_id": sample_id, "read_type": read_type})
        rows[(sample_id, read_type)].update(
            {
                "clean_reads": to_int(record.get("num_seqs"), ""),
                "clean_bases": to_int(record.get("sum_len"), ""),
                "clean_mean_length": to_float(record.get("avg_len"), ""),
                "clean_mean_quality": to_float(record.get("AvgQual"), ""),
            }
        )

    table_data = {}
    for key in sorted(rows):
        row = rows[key]
        table_data[row_label(row["sample_id"], row["read_type"])] = {
            "sample_id": row["sample_id"],
            "read_type": row["read_type"],
            "raw_reads": row.get("raw_reads", ""),
            "raw_bases": row.get("raw_bases", ""),
            "raw_mean_length": row.get("raw_mean_length", ""),
            "raw_mean_quality": row.get("raw_mean_quality", ""),
            "clean_reads": row.get("clean_reads", ""),
            "clean_bases": row.get("clean_bases", ""),
            "clean_mean_length": row.get("clean_mean_length", ""),
            "clean_mean_quality": row.get("clean_mean_quality", ""),
        }

    return {
        "id": "h2seq_read_statistics",
        "section_name": "Read Statistics",
        "description": "Raw and quality-controlled read statistics harmonised across long- and short-read inputs.",
        "plot_type": "table",
        "pconfig": {
            "id": "h2seq_read_statistics_table",
            "title": "Read Statistics",
        },
        "headers": {
            "sample_id": {"title": "Sample ID"},
            "read_type": {"title": "Read Type"},
            "raw_reads": {"title": "Raw reads"},
            "raw_bases": {"title": "Raw bases"},
            "raw_mean_length": {"title": "Raw mean length", "format": "{:,.1f}"},
            "raw_mean_quality": {"title": "Raw mean Q", "format": "{:,.2f}"},
            "clean_reads": {"title": "Reads after QC"},
            "clean_bases": {"title": "Bases after QC"},
            "clean_mean_length": {"title": "Clean mean length", "format": "{:,.1f}"},
            "clean_mean_quality": {"title": "Clean mean Q", "format": "{:,.2f}"},
        },
        "data": table_data,
    }


def build_variant_section(outdir):
    rows = {}

    for path in sorted(outdir.rglob("*.coverage_summary.tsv")):
        for record in load_tsv_rows(path):
            key = (record["sample_id"], record["read_type"])
            rows.setdefault(
                key,
                {
                    "sample_id": record["sample_id"],
                    "read_type": record["read_type"],
                    "reference": record.get("reference_name", ""),
                    "snp_number": 0,
                    "indel_number": 0,
                },
            )

    for path in sorted(outdir.rglob("*.annotated.tsv")):
        for record in load_tsv_rows(path):
            key = (record["sample_id"], record["read_type"])
            row = rows.setdefault(
                key,
                {
                    "sample_id": record["sample_id"],
                    "read_type": record["read_type"],
                    "reference": record.get("chrom", ""),
                    "snp_number": 0,
                    "indel_number": 0,
                },
            )
            if not row["reference"]:
                row["reference"] = record.get("chrom", "")
            variant_type = (record.get("type") or "").upper()
            if variant_type == "SNP":
                row["snp_number"] += 1
            elif variant_type == "INDEL":
                row["indel_number"] += 1

    table_data = {}
    for key in sorted(rows):
        row = rows[key]
        table_data[row_label(row["sample_id"], row["read_type"])] = {
            "sample_id": row["sample_id"],
            "read_type": row["read_type"],
            "reference": row["reference"],
            "snp_number": row["snp_number"],
            "indel_number": row["indel_number"],
        }

    return {
        "id": "h2seq_variant_calling",
        "section_name": "Variant Calling",
        "description": "Counts of filtered variants contributing to the final consensus sequences.",
        "plot_type": "table",
        "pconfig": {
            "id": "h2seq_variant_calling_table",
            "title": "Variant Calling",
        },
        "headers": {
            "sample_id": {"title": "Sample ID"},
            "read_type": {"title": "Read Type"},
            "reference": {"title": "Reference"},
            "snp_number": {"title": "SNP number"},
            "indel_number": {"title": "Indel number"},
        },
        "data": table_data,
    }


def write_json(path, payload):
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    write_json(args.coverage_output, build_coverage_section(outdir))
    write_json(args.read_output, build_read_section(outdir))
    write_json(args.variant_output, build_variant_section(outdir))


if __name__ == "__main__":
    main()
