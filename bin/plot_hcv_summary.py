#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path


import matplotlib.pyplot as plt


FEATURES = [
    "Genome",
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
    parser = argparse.ArgumentParser(description="Create HCV summary plots.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read-type", required=True, choices=["long", "short"])
    parser.add_argument("--reference-name", required=True)
    parser.add_argument("--per-base-bed-gz", required=True)
    parser.add_argument("--coverage-summary", required=True)
    parser.add_argument("--hcv-coverage", required=True)
    parser.add_argument("--min-depth", required=True, type=float)
    parser.add_argument("--depth-plot-output", required=True)
    parser.add_argument("--feature-plot-output", required=True)
    return parser.parse_args()


def load_coverage_summary(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if len(rows) != 1:
        raise ValueError(f"Expected one row in {path}, found {len(rows)}")
    row = rows[0]
    return {
        "reference_length": int(float(row["reference_length"])),
        "genome_coverage_pct": float(row["genome_coverage_pct"]),
        "mean_depth": float(row["mean_depth"]),
    }


def load_per_base_segments(path):
    segments = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            chrom, start, end, depth = line.split("\t")[:4]
            segments.append((int(start), int(end), float(depth)))
    return segments


def load_hcv_feature_coverage(path):
    values = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw = (row.get("coverage_pct") or "").strip()
            values[row["feature"]] = float(raw) if raw else 0.0
    return values


def build_depth_points(segments, reference_length):
    points_x = [0]
    points_y = [0.0]
    last_end = 0
    for start, end, depth in segments:
        if start > last_end:
            points_x.extend([start, start])
            points_y.extend([0.0, 0.0])
        points_x.extend([start, end])
        points_y.extend([depth, depth])
        last_end = end
    if last_end < reference_length:
        points_x.append(reference_length)
        points_y.append(0.0)
    return points_x, points_y


def create_depth_plot(args, coverage_summary, segments):
    reference_length = coverage_summary["reference_length"]
    x_values, y_values = build_depth_points(segments, reference_length)
    true_max_depth = max(y_values) if y_values else 0.0
    y_max = min(max(true_max_depth, args.min_depth, 1.0), 1000.0)

    fig, ax = plt.subplots(figsize=(10.5, 3.6), dpi=300)
    ax.plot(x_values, [min(value, y_max) for value in y_values], color="#127a5a", linewidth=1.2)
    ax.axhline(args.min_depth, color="#d9472b", linestyle="--", linewidth=1.0)
    ax.text(reference_length * 0.995, args.min_depth + (y_max * 0.015), f"{args.min_depth:g}x threshold", ha="right", va="bottom", fontsize=8)
    ax.set_title(f"{args.sample_id} ({args.read_type}) coverage across selected reference genome", fontsize=11, weight="bold")
    ax.set_ylabel("Depth", labelpad=16)
    ax.set_xlabel("Reference position")
    ax.set_xlim(0, reference_length)
    tick_positions = [pos for pos in range(0, 10001, 1000) if pos <= reference_length]
    ax.set_xticks(tick_positions)
    ax.set_ylim(0, y_max)
    ax.tick_params(axis="y", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    if true_max_depth > y_max:
        fig.text(
            0.01,
            0.01,
            f"Note: maximum depth capped at {int(y_max)}x for plotting rather than the true maximum of {true_max_depth:.0f}x",
            ha="left",
            va="bottom",
            fontsize=7,
        )
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(args.depth_plot_output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def create_feature_plot(args, coverage_summary, feature_coverage):
    values = {"Genome": coverage_summary["genome_coverage_pct"]}
    values.update(feature_coverage)
    heights = [values.get(feature, 0.0) for feature in FEATURES]

    fig, ax = plt.subplots(figsize=(10.5, 3.6), dpi=300)
    bars = ax.bar(range(len(FEATURES)), heights, color="#5d84b8")
    ax.set_ylim(0, 100)
    ax.set_ylabel("Coverage (%)")
    ax.set_xticks(range(len(FEATURES)))
    ax.set_xticklabels(FEATURES, rotation=45, ha="right")
    ax.set_title("Genome and HCV feature coverage", fontsize=11, weight="bold", pad=16)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    for bar, height in zip(bars, heights):
        ax.text(bar.get_x() + bar.get_width() / 2.0, min(height + 1.5, 99.5), f"{height:.1f}", ha="center", va="bottom", fontsize=6)
    fig.tight_layout()
    fig.savefig(args.feature_plot_output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    coverage_summary = load_coverage_summary(args.coverage_summary)
    segments = load_per_base_segments(args.per_base_bed_gz)
    feature_coverage = load_hcv_feature_coverage(args.hcv_coverage)
    create_depth_plot(args, coverage_summary, segments)
    create_feature_plot(args, coverage_summary, feature_coverage)


if __name__ == "__main__":
    main()
