#!/usr/bin/env python
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Select the best reference genome based on kallisto-derived abundances.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input TSV file')
    parser.add_argument('-s', '--sample_name', type=str, required=True, help='Sample name')
    parser.add_argument('--top', type=float, required=True, help='Top percent to select target_ids within')
    parser.add_argument('-b', '--best_ref_txt', required=True, type=str, default='best_reference.txt', help='Output TXT file path (will have just the reference ID)')
    parser.add_argument('--salmon', action='store_true', help='Salmon was used to generate the input (rather than kallisto)')
    parser.add_argument('-o', '--output', type=str, default='best_reference.tsv', help='Output TSV file path')
    return parser.parse_args()

def process_data(input_path, salmon, sample_name, top_percent, best_ref_txt, output_path):
    df = pd.read_csv(input_path, sep='\t')
    
    # update headers based on tool used
    if salmon:
        new_headers = ['target_id', 'length', 'eff_length', 'tpm','est_counts']
        df.columns = new_headers
    
    # Choose the best target_id based on the highest tpm
    best_row = df.loc[df['tpm'].idxmax()]
    best_tpm = best_row['tpm']
    best_target_id = best_row['target_id']
    
    # write the best target id
    with open(best_ref_txt, 'w') as f:
        f.write(f'{best_target_id}\n')
    
    # Extract genotype and subtype
    best_genotype, best_subtype = extract_genotype_subtype(best_target_id)
    
    # Find the n-best target_ids within top_percent of the highest tpm
    threshold_tpm = best_tpm - (best_tpm * (top_percent / 100))
    close_hits_df = df[df['tpm'] >= threshold_tpm]
    close_hits = ';'.join(close_hits_df['target_id'].tolist())
    
    # Write out the TSV file
    with open(output_path, 'w') as f:
        f.write(f'sample_id\tgenotype\tsubtype\tbest_ref\tclose_hits\n')
        f.write(f'{sample_name}\t{best_genotype}\t{best_subtype}\t{best_target_id}\t{close_hits}\n')

def extract_genotype_subtype(target_id):
    parts = target_id.split('_')
    genotype = ''.join([char for char in parts[0] if char.isdigit()])
    subtype = parts[0][len(genotype):]
    return genotype, subtype

if __name__ == "__main__":
    args = parse_args()
    process_data(args.input, args.salmon, args.sample_name, args.top, args.best_ref_txt, args.output)
