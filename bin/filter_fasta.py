#!/usr/bin/env python

import sys
from Bio import SeqIO
import os

def output_stem(prefix, name):
    return name if name.startswith(f"{prefix}.") or name == prefix else f"{prefix}.{name}"

def filter_sequences(input_fastas, prefix):
    for f in input_fastas:
        name = os.path.splitext(os.path.basename(f))[0]
        stem = output_stem(prefix, name)

        reporting_file = f"{stem}.reporting.fasta"
        empty_file = f"{stem}.empty.fasta"

        with open(reporting_file, 'w') as reporting, open(empty_file, 'w') as empty:
            for record in SeqIO.parse(f, "fasta"):
                if all(base == 'N' for base in record.seq):
                    SeqIO.write(record, empty, "fasta")
                else:
                    SeqIO.write(record, reporting, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: filter_fasta.py prefix file.fasta [file.fasta ...]")
        sys.exit(1)

    prefix = sys.argv[1]
    input_fastas = sys.argv[2:]

    filter_sequences(input_fastas, prefix)
