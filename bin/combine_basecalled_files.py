#!/usr/bin/env python

import argparse
import csv
import gzip
import os
from pathlib import Path
import sys

def validate_input(row, base_dir):
    required_fields = ['id', 'run_name', 'barcode']
    for field in required_fields:
        if field not in row or not row[field]:
            raise ValueError(f"{field} is mandatory and missing in row: {row}")
    location = list(base_dir.glob(f"**/fastq_pass/{row['barcode']}"))
    if not location:
        raise ValueError(f"Barcode path does not exist: {row['barcode']} in {base_dir}")
    if not os.path.isdir(location[0]):
        raise ValueError(f"Barcode path is not a directory: {location[0]}")

def combine_fastq_files(fastq_files, output_file, overwrite, compressed):
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    if not overwrite and Path(output_file).exists():
        error_message = f"[ \033[91mERROR\033[0m ] : FileExistsError: Output file {output_file} already exists. \n\nTip: Use option '--overwrite' when running the script to enable overwriting of existing files."
        print(error_message, file=sys.stderr)
        sys.exit(1)

    if compressed:
        with gzip.open(output_file, 'wt') as outfile:
            for fname in fastq_files:
                with gzip.open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
    else:
        with open(output_file, 'w') as outfile:
            for fname in fastq_files:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

def read_csv(filepath):
    with open(filepath, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
    return rows

def update_barcode(barcode):
    if barcode.startswith('BC'):
        return 'barcode' + barcode[2:]
    return barcode

def parse_final_summary(summary_file):
    instrument, flow_cell_id, protocol = None, None, None
    with open(summary_file) as file:
        for line in file:
            if line.startswith("instrument="):
                instrument = line.split("=")[1].strip()
            elif line.startswith("flow_cell_id="):
                flow_cell_id = line.split("=")[1].strip()
            elif line.startswith("protocol="):
                protocol = line.split("=")[1].strip()
    return instrument, flow_cell_id, protocol

def main(basecall_default_path, input_path, legacy, overwrite, outdir, compressed, csv_name):
    # Read the spreadsheet
    rows = read_csv(input_path)
    csv_rows = []

    for row in rows:
        if not legacy:
            if 'barcode' not in row:
                raise KeyError(f"'barcode' column is missing in row: {row}")
            row['barcode'] = update_barcode(row['barcode'])
        
        id = row['id']
        run_name = row['run_name']
        barcode_path = row['barcode']

        base_dir = Path(basecall_default_path) / run_name
        validate_input(row, base_dir)

        barcode_dir = list(base_dir.glob(f"**/fastq_pass/{Path(barcode_path).name}"))
        if not barcode_dir:
            raise ValueError(f"Barcode directory not found for {barcode_path} in {base_dir}")

        parent_dir = barcode_dir[0].parent.parent
        summary_file = list(parent_dir.glob('final_summary_*.txt'))
        if not summary_file:
            raise ValueError(f"No final_summary_*.txt file found in {parent_dir}")
        
        instrument, flow_cell_id, protocol = parse_final_summary(summary_file[0])

        fastq_dir = barcode_dir[0]
        fastq_files = [f for f in fastq_dir.glob('*.fastq') if not f.name.startswith('.')]
        if not fastq_files:
            raise ValueError(f"No .fastq files found in {fastq_dir}")

        long_reads = Path(outdir) / f"{id}.fastq.gz" if compressed else Path(outdir) / f"{id}.fastq"
        combine_fastq_files(fastq_files, long_reads, overwrite, compressed)
        
        csv_rows.append({
            'sample': id,
            'long_reads': long_reads,
            'short_reads_1': '',
            'short_reads_2': '',
            'run_name': run_name,
            'barcode': barcode_path,
            'instrument': instrument,
            'flow_cell_id': flow_cell_id,
            'protocol': protocol
        })

        print(f"[ Sample: \033[92m{id}\033[0m ] : combined {len(fastq_files)} fastq files into {long_reads}")

    with open(csv_name, 'w', newline='') as csvfile:
        fieldnames = ['sample', 'long_reads', 'short_reads_1', 'short_reads_2', 'run_name', 'barcode', 'instrument', 'flow_cell_id', 'protocol']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine all output .fastq files from a guppy/dorado basecalling run into a single file per id",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", help="Input spreadsheet path")
    parser.add_argument("-p", "--basecall_default_path", default="/var/lib/minknow/data", help="Path to basecalling base directory")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory for combined files")
    parser.add_argument("-c", "--csv_name", default="h2seq_pipeline_input.csv", help="Name of the output CSV file")
    parser.add_argument("--legacy", action="store_true", help="Use legacy barcode names")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--compressed", action="store_true", help="Input fastq files are compressed (output will also be compressed)")
    
    args = parser.parse_args()
    main(args.basecall_default_path, args.input, args.legacy, args.overwrite, args.outdir, args.compressed, args.csv_name)
