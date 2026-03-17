#!/usr/bin/env python
import argparse
from statistics import mean, stdev
import sys

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate mean read length and standard deviation for a FASTQ file."
    )
    parser.add_argument("fastq_file", help="Input FASTQ file")
    return parser.parse_args()


def calculate_read_length_stats(fastq_file):
    read_lengths = [len(record.seq) for record in SeqIO.parse(fastq_file, "fastq")]

    if not read_lengths:
        raise ValueError(f"No reads were found in FASTQ file: {fastq_file}")

    mean_length = round(mean(read_lengths))
    std_dev = round(stdev(read_lengths)) if len(read_lengths) > 1 else 0
    return mean_length, std_dev


def write_value(value, filename):
    with open(filename, "w") as file:
        file.write(f"{value}\n")


if __name__ == "__main__":
    args = parse_args()
    try:
        mean_length, std_dev = calculate_read_length_stats(args.fastq_file)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)

    write_value(mean_length, "mean_length.txt")
    write_value(std_dev, "std_dev.txt")
