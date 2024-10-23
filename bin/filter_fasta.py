#!/usr/bin/env python

import sys
from Bio import SeqIO
import os

def filter_sequences(input_fastas, prefix):
    for f in input_fastas:
        name = os.path.splitext(os.path.basename(f))[0]

        reporting_file = f"{name}.reporting.fasta"
        empty_file = f"{name}.empty.fasta"

        with open(reporting_file, 'w') as reporting, open(empty_file, 'w') as empty:
            for record in SeqIO.parse(f, "fasta"):
                if all(base == 'N' for base in record.seq):
                    SeqIO.write(record, empty, "fasta")
                else:
                    SeqIO.write(record, reporting, "fasta")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: filter_sequences.py name file.fasta")
        sys.exit(1)

    prefix = sys.argv[1]
    input_fastas = sys.argv[2:]

    filter_sequences(input_fastas, prefix)
