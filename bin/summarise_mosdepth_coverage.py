#!/usr/bin/env python3
import argparse
import csv
import gzip


def parse_args():
    parser = argparse.ArgumentParser(description="Summarise whole-genome coverage from mosdepth output.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read-type", required=True, choices=["long", "short"])
    parser.add_argument("--reference-name", required=True)
    parser.add_argument("--genome-bed-gz", required=True)
    parser.add_argument("--per-base-bed-gz", required=True)
    parser.add_argument("--min-depth", required=True, type=float)
    parser.add_argument("--summary-output", required=True)
    return parser.parse_args()


def read_mosdepth_regions(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                raise ValueError(f"Unexpected mosdepth regions format in {path}: {line}")
            yield {
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "name": fields[3],
                "mean_depth": float(fields[4]),
            }


def read_per_base_segments(path):
    segments = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"Unexpected mosdepth per-base format in {path}: {line}")
            segments.append((fields[0], int(fields[1]), int(fields[2]), float(fields[3])))
    return segments


def interval_coverage(interval, segments, min_depth):
    total = max(interval["end"] - interval["start"], 0)
    covered = 0
    depth_sum = 0.0

    for chrom, seg_start, seg_end, depth in segments:
        if chrom != interval["chrom"]:
            continue
        start = max(interval["start"], seg_start)
        end = min(interval["end"], seg_end)
        if end <= start:
            continue
        length = end - start
        if depth >= min_depth:
            covered += length
        depth_sum += length * depth

    coverage_pct = (covered / total * 100.0) if total else 0.0
    mean_depth = (depth_sum / total) if total else 0.0
    return total, covered, coverage_pct, mean_depth


def main():
    args = parse_args()

    genome_rows = list(read_mosdepth_regions(args.genome_bed_gz))
    if not genome_rows:
        raise ValueError(f"No whole-genome mosdepth rows found in {args.genome_bed_gz}")
    per_base_segments = read_per_base_segments(args.per_base_bed_gz)

    with open(args.summary_output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_id",
                "read_type",
                "reference_name",
                "reference_length",
                "positions_covered",
                "genome_coverage_pct",
                "mean_depth",
            ]
        )
        for genome_row in genome_rows:
            genome_interval = {
                "chrom": genome_row["chrom"],
                "start": genome_row["start"],
                "end": genome_row["end"],
                "name": genome_row["name"],
            }
            reference_length, positions_covered, genome_coverage_pct, mean_depth = interval_coverage(
                genome_interval, per_base_segments, args.min_depth
            )
            writer.writerow(
                [
                    args.sample_id,
                    args.read_type,
                    genome_row["chrom"],
                    reference_length,
                    positions_covered,
                    genome_coverage_pct,
                    mean_depth,
                ]
            )


if __name__ == "__main__":
    main()
