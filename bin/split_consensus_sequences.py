#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys

def extract_sequences(pattern_file, fasta_file, prefix):
    # Read the pattern
    try:
        with open(pattern_file) as f:
            pattern = f.readline().strip()
    except FileNotFoundError:
        sys.exit(f"Error: Pattern file '{pattern_file}' not found.")
    except Exception as e:
        sys.exit(f"Error reading pattern file: {e}")

    # Initialize storage
    main_sequences = []
    alt_sequences = []

    # Parse FASTA file
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            sys.exit(f"Error: No valid sequences found in '{fasta_file}'")
    except FileNotFoundError:
        sys.exit(f"Error: FASTA file '{fasta_file}' not found.")
    except Exception as e:
        sys.exit(f"Error parsing FASTA file: {e}")

    # Sort sequences into main and alternative groups
    for record in records:
        if pattern in record.description:
            main_sequences.append(record)
        else:
            alt_sequences.append(record)

    # Write main sequences if any exist
    if main_sequences:
        main_output = f"{prefix}.consensus_main.fa"
        try:
            SeqIO.write(main_sequences, main_output, "fasta")
            print(f"Wrote {len(main_sequences)} sequence(s) to {main_output}")
        except Exception as e:
            sys.exit(f"Error writing main sequences: {e}")
    else:
        print("No sequences matched the pattern - no main reference file created")

    # Write alternative sequences only if they exist
    if alt_sequences:
        for i, record in enumerate(alt_sequences, 1):
            alt_output = f"{prefix}.consensus_alt{i}.fa"
            try:
                SeqIO.write([record], alt_output, "fasta")
                print(f"Wrote sequence to {alt_output}")
            except Exception as e:
                sys.exit(f"Error writing alternative sequence {i}: {e}")
    else:
        print("No alternative sequences found - no alt reference files created")

def main():
    parser = argparse.ArgumentParser(description='Extract sequences from a FASTA file based on a pattern.')
    parser.add_argument('pattern_file', help='File containing the pattern to match')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--prefix', default='custom_id', help='Prefix for output files (default: custom_id)')

    args = parser.parse_args()
    extract_sequences(args.pattern_file, args.fasta_file, args.prefix)

if __name__ == '__main__':
    main()