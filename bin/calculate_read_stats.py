#!/usr/bin/env python
from Bio import SeqIO
from statistics import mean, stdev
import sys

def calculate_read_length_stats(fastq_file):
    read_lengths = []
    for record in SeqIO.parse(fastq_file, "fastq"):
        read_lengths.append(len(record.seq))
    mean_length = round(mean(read_lengths))
    std_dev = round(stdev(read_lengths))
    return mean_length, std_dev

def write_value(value, filename):
    with open(filename, 'w') as file:
        file.write(str(value) + '\n')

if __name__ == "__main__":
    fastq_file = sys.argv[1]
    mean_length, std_dev = calculate_read_length_stats(fastq_file)
    write_value(mean_length, "mean_length.txt")
    write_value(std_dev, "std_dev.txt")
