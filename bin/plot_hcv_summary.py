#!/usr/bin/env python3
import argparse
import csv
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
    parser.add_argument("--coverage-summary", required=True)
    parser.add_argument("--hcv-coverage", required=True)
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
        "genome_coverage_pct": float(row["genome_coverage_pct"]),
    }


def load_hcv_feature_coverage(path):
    values = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw = (row.get("coverage_pct") or "").strip()
            values[row["feature"]] = float(raw) if raw else 0.0
    return values


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
    feature_coverage = load_hcv_feature_coverage(args.hcv_coverage)
    create_feature_plot(args, coverage_summary, feature_coverage)


if __name__ == "__main__":
    main()
