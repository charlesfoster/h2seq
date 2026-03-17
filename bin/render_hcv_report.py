#!/usr/bin/env python3
import argparse
import csv
import textwrap
from datetime import date

import matplotlib.image as mpimg
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Render a one-page HCV PDF report.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read-type", required=True, choices=["long", "short"])
    parser.add_argument("--best-reference-tsv", required=True)
    parser.add_argument("--coverage-summary", required=True)
    parser.add_argument("--depth-plot", required=True)
    parser.add_argument("--feature-plot", required=True)
    parser.add_argument("--logo", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--reference-selection-tool", required=True)
    parser.add_argument("--consensus-min-depth", required=True)
    parser.add_argument("--snv-min-af", required=True)
    parser.add_argument("--indel-min-af", required=True)
    parser.add_argument("--long-reads-min-len", required=True)
    parser.add_argument("--long-reads-max-len", required=True)
    parser.add_argument("--short-reads-min-len", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_single_tsv_row(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if len(rows) != 1:
        raise ValueError(f"Expected one row in {path}, found {len(rows)}")
    return rows[0]


def build_summary_sentence(args):
    snv_af_pct = float(args.snv_min_af) * 100.0
    indel_af_pct = float(args.indel_min_af) * 100.0
    if args.read_type == "long":
        qc_description = (
            f"Raw reads were quality controlled to retain sequences within the configured long-read length range "
            f"of {args.long_reads_min_len} to {args.long_reads_max_len} bases and to exclude lower-quality reads."
        )
    else:
        qc_description = (
            f"Raw reads were quality controlled to retain reads meeting the configured short-read filtering thresholds, "
            f"including a minimum read length of {args.short_reads_min_len} bases."
        )
    version = args.pipeline_version or "unknown"
    return (
        f"{args.sample_id} was analysed on {date.today().isoformat()} using version v{version} of the H2seq bioinformatics pipeline. "
        f"{qc_description} A closest reference was selected from the configured reference panel based on a read mapping approach, "
        f"reads were aligned to that reference, SNVs were retained from a minimum depth of {args.consensus_min_depth} and minimum allele frequency of {snv_af_pct:.1f}%, "
        f"indels were retained using a minimum allele frequency of {indel_af_pct:.1f}%, and a consensus genome was assembled using a minimum consensus depth of {args.consensus_min_depth}. "
        f"Coding-region coverage metrics were taken from a downstream resistance-analysis report."
    )


def main():
    args = parse_args()
    best_ref = read_single_tsv_row(args.best_reference_tsv)
    coverage = read_single_tsv_row(args.coverage_summary)
    logo = mpimg.imread(args.logo)
    depth_plot = mpimg.imread(args.depth_plot)
    feature_plot = mpimg.imread(args.feature_plot)
    summary_sentence = build_summary_sentence(args)

    fig = plt.figure(figsize=(8.27, 11.69), dpi=200)
    gs = fig.add_gridspec(nrows=28, ncols=12, left=0.05, right=0.95, top=0.97, bottom=0.04, hspace=0.45, wspace=0.4)

    ax_title = fig.add_subplot(gs[0:2, 0:8])
    ax_title.axis("off")
    ax_title.text(0.0, 0.7, "H2seq Workflow Results", fontsize=18, fontweight="bold", ha="left", va="center")

    ax_logo = fig.add_subplot(gs[0:2, 9:12])
    ax_logo.imshow(logo)
    ax_logo.axis("off")

    ax_table = fig.add_subplot(gs[2:7, 0:12])
    ax_table.axis("off")
    table_rows = [
        ["Date", date.today().isoformat()],
        ["Sample", args.sample_id],
        ["Genotype", best_ref.get("genotype", "")],
        ["Subtype", best_ref.get("subtype", "")],
        ["Closest reference", best_ref.get("best_ref", "")],
        ["Genome coverage", f"{float(coverage.get('genome_coverage_pct', 0.0)):.1f}%"],
        ["Mean depth", f"{float(coverage.get('mean_depth', 0.0)):.1f}x"],
    ]
    table = ax_table.table(cellText=table_rows, colWidths=[0.25, 0.75], cellLoc="left", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#aeb8c4")
        if col == 0:
            cell.set_facecolor("#dfe8f2")
            cell.set_text_props(weight="bold")
        else:
            cell.set_facecolor("#f7f9fc")

    ax_plot1 = fig.add_subplot(gs[7:14, 0:12])
    ax_plot1.imshow(depth_plot)
    ax_plot1.axis("off")

    ax_plot2 = fig.add_subplot(gs[14:21, 0:12])
    ax_plot2.imshow(feature_plot)
    ax_plot2.axis("off")

    ax_footer = fig.add_subplot(gs[21:28, 0:12])
    ax_footer.axis("off")
    ax_footer.text(0.0, 0.98, "Brief Description", fontsize=11, fontweight="bold", ha="left", va="top")
    ax_footer.text(0.0, 0.86, textwrap.fill(summary_sentence, width=120), fontsize=9, ha="left", va="top")

    fig.savefig(args.output, format="pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
