#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt

MOSDEPTH_COVERAGE_THRESHOLDS = [1, 5, 10, 20, 30, 50, 100]
HCV_FEATURES = [
    "Polyprotein",
    "Core",
    "E1",
    "E2",
    "p7",
    "NS2",
    "NS3",
    "NS4A",
    "NS4B",
    "NS5A",
    "NS5B",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Build custom MultiQC sections from pipeline outputs.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--coverage-output", required=True)
    parser.add_argument("--region-coverage-output", required=True)
    parser.add_argument("--region-coverage-pdf-output", required=True)
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


def summarise_seqkit_rows(path):
    rows = load_tsv_rows(path)
    if not rows:
        raise ValueError(f"No seqkit rows found in {path}")

    counts = [to_int(row.get("num_seqs"), 0) for row in rows]
    bases = [to_int(row.get("sum_len"), 0) for row in rows]
    weighted_quality_num = 0.0
    weighted_quality_den = 0

    for row, base_count in zip(rows, bases):
        avg_qual = to_float(row.get("AvgQual"))
        if avg_qual is not None and base_count:
            weighted_quality_num += avg_qual * base_count
            weighted_quality_den += base_count

    if len(set(counts)) == 1:
        read_count = counts[0]
    else:
        read_count = sum(counts)

    total_bases = sum(bases)
    mean_length = (total_bases / sum(counts)) if sum(counts) else None
    mean_quality = (weighted_quality_num / weighted_quality_den) if weighted_quality_den else None

    return {
        "num_seqs": read_count,
        "sum_len": total_bases,
        "avg_len": mean_length,
        "AvgQual": mean_quality,
    }


def to_float(value, default=None):
    if value in {"", None}:
        return default
    return float(value)


def to_int(value, default=None):
    if value in {"", None}:
        return default
    return int(float(value))


def parse_mosdepth_global_dist(path):
    coverage_pct_by_threshold = {depth: None for depth in MOSDEPTH_COVERAGE_THRESHOLDS}
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
            if depth in coverage_pct_by_threshold:
                coverage_pct_by_threshold[depth] = fraction * 100.0
            if fraction >= 0.5 and depth > median_depth:
                median_depth = depth

    summary = {"median_coverage": median_depth}
    for depth, value in coverage_pct_by_threshold.items():
        summary[f"coverage_{depth}x_pct"] = value if value is not None else 0.0
    return summary


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
            "median_coverage": row.get("median_coverage", ""),
            "mean_coverage": row.get("mean_coverage", ""),
            "min_coverage": row.get("min_coverage", ""),
            "max_coverage": row.get("max_coverage", ""),
            "genome_length": row.get("genome_length", ""),
        }
        for depth in MOSDEPTH_COVERAGE_THRESHOLDS:
            table_data[row_label(row["sample_id"], row["read_type"])][f"coverage_{depth}x_pct"] = row.get(
                f"coverage_{depth}x_pct", ""
            )

    headers = {}
    for depth in MOSDEPTH_COVERAGE_THRESHOLDS:
        headers[f"coverage_{depth}x_pct"] = {
            "title": f"≥ {depth}X",
            "format": "{:,.1f}",
            "suffix": "%",
            "hidden": depth != 10,
        }
    headers.update(
        {
            "median_coverage": {"title": "Median", "format": "{:,.1f}X"},
            "mean_coverage": {"title": "Mean Cov.", "format": "{:,.1f}"},
            "min_coverage": {"title": "Min Cov.", "format": "{:,.1f}"},
            "max_coverage": {"title": "Max Cov.", "format": "{:,.1f}"},
            "genome_length": {"title": "Genome length"},
        }
    )

    return {
        "id": "h2seq_coverage_statistics",
        "section_name": "Coverage Statistics",
        "description": "Coverage summary metrics derived from mosdepth outputs.",
        "plot_type": "table",
        "pconfig": {
            "id": "h2seq_coverage_statistics_table",
            "title": "Coverage Statistics",
        },
        "headers": headers,
        "data": table_data,
    }


def build_read_section(outdir):
    rows = {}

    for path in sorted(outdir.rglob("*.raw_long.tsv")) + sorted(outdir.rglob("*.raw_short.tsv")):
        sample_id, read_type = infer_sample_and_read_type(path)
        record = summarise_seqkit_rows(path)
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
        record = summarise_seqkit_rows(path)
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


def build_region_coverage_rows(outdir):
    rows = {}

    combined_summary = outdir / "combined_results_summary.csv"
    if combined_summary.exists():
        with open(combined_summary, newline="") as handle:
            for record in csv.DictReader(handle):
                sample_id = record.get("sample_id", "")
                read_type = record.get("read_type", "")
                if not sample_id or not read_type:
                    continue
                key = (sample_id, read_type)
                rows.setdefault(
                    key,
                    {
                        "sample_id": sample_id,
                        "read_type": read_type,
                        **{feature: 0.0 for feature in HCV_FEATURES},
                    },
                )

    if not rows:
        for path in sorted(outdir.rglob("*.coverage_summary.tsv")):
            for record in load_tsv_rows(path):
                sample_id = record.get("sample_id", "")
                read_type = record.get("read_type", "")
                if not sample_id or not read_type:
                    continue
                key = (sample_id, read_type)
                rows.setdefault(
                    key,
                    {
                        "sample_id": sample_id,
                        "read_type": read_type,
                        **{feature: 0.0 for feature in HCV_FEATURES},
                    },
                )

    for path in sorted(outdir.rglob("*.hcv_glue_coverage.tsv")):
        for record in load_tsv_rows(path):
            sample_id = record.get("sample_id", "")
            read_type = record.get("read_type", "")
            feature = record.get("feature", "")
            if feature not in HCV_FEATURES:
                continue
            key = (sample_id, read_type)
            rows.setdefault(
                key,
                {
                    "sample_id": sample_id,
                    "read_type": read_type,
                    **{feature_name: 0.0 for feature_name in HCV_FEATURES},
                },
            )
            rows[key][feature] = to_float(record.get("coverage_pct"), 0.0)

    return rows


def build_region_coverage_section(outdir):
    rows = build_region_coverage_rows(outdir)

    if not rows:
        return {}

    ordered_rows = [rows[key] for key in sorted(rows)]
    sample_labels = [row_label(row["sample_id"], row["read_type"]) for row in ordered_rows]
    matrix = [[row.get(feature, 0.0) for feature in HCV_FEATURES] for row in ordered_rows]

    return {
        "id": "h2seq_region_coverage",
        "section_name": "Coverage Per Genomic Region",
        "description": "HCV gene and polyprotein coverage percentages parsed from HCV-GLUE reports.",
        "plot_type": "heatmap",
        "pconfig": {
            "id": "h2seq_region_coverage_heatmap",
            "title": "Coverage Per Genomic Region",
            "xlab": "Genomic region",
            "ylab": "Sample",
            "zlab": "Coverage (%)",
            "min": 0,
            "max": 100,
            "xcats_samples": False,
            "ycats_samples": True,
        },
        "data": matrix,
        "xcats": HCV_FEATURES,
        "ycats": sample_labels,
    }


def write_region_coverage_heatmap(path, outdir):
    rows = build_region_coverage_rows(outdir)

    if not rows:
        fig, ax = plt.subplots(figsize=(6, 1.8), dpi=300)
        ax.axis("off")
        ax.text(0.5, 0.5, "No HCV-GLUE region coverage available", ha="center", va="center", fontsize=11)
        fig.tight_layout()
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        return

    ordered_rows = [rows[key] for key in sorted(rows)]
    sample_labels = [row_label(row["sample_id"], row["read_type"]) for row in ordered_rows]
    matrix = [[row.get(feature, 0.0) if row.get(feature, None) is not None else float("nan") for feature in HCV_FEATURES] for row in ordered_rows]

    fig_height = max(2.6, 0.55 * len(sample_labels) + 1.8)
    fig, ax = plt.subplots(figsize=(10.5, fig_height), dpi=300)
    cmap = plt.cm.YlGnBu.copy()
    cmap.set_bad(color="#f2f2f2")
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=100)
    ax.set_title("Coverage Per Genomic Region", fontsize=13, weight="bold", pad=12)
    ax.set_xticks(range(len(HCV_FEATURES)))
    ax.set_xticklabels(HCV_FEATURES, rotation=45, ha="right")
    ax.set_yticks(range(len(sample_labels)))
    ax.set_yticklabels(sample_labels)

    for row_idx, row in enumerate(matrix):
        for col_idx, value in enumerate(row):
            if value == value:
                ax.text(col_idx, row_idx, f"{value:.1f}", ha="center", va="center", fontsize=7, color="black")

    cbar = fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Coverage (%)")
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_json(path, payload):
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    write_json(args.coverage_output, build_coverage_section(outdir))
    write_json(args.region_coverage_output, build_region_coverage_section(outdir))
    write_region_coverage_heatmap(args.region_coverage_pdf_output, outdir)
    write_json(args.read_output, build_read_section(outdir))
    write_json(args.variant_output, build_variant_section(outdir))


if __name__ == "__main__":
    main()
