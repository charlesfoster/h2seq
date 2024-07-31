#!/usr/bin/env python

import sys
from Bio import SeqIO

def filter_sequences(input_fasta, prefix):
    reporting_file = f"{prefix}.reporting.fasta"
    empty_file = f"{prefix}.empty.fasta"

    with open(reporting_file, 'w') as reporting, open(empty_file, 'w') as empty:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if all(base == 'N' for base in record.seq):
                SeqIO.write(record, empty, "fasta")
            else:
                SeqIO.write(record, reporting, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: filter_sequences.py file.fasta name")
        sys.exit(1)

    input_fasta = sys.argv[1]
    prefix = sys.argv[2]

    filter_sequences(input_fasta, prefix)
