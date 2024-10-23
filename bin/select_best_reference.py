#!/usr/bin/env python
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Select the best reference genome based on kallisto-derived abundances.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input TSV file')
    parser.add_argument('-s', '--sample_name', type=str, required=True, help='Sample name')
    parser.add_argument('--top', type=float, required=True, help='Top percent of TPM to select target_ids within')
    parser.add_argument('-b', '--best_ref_txt', required=True, type=str, default='best_reference.txt', help='Output TXT file path (will have just the reference ID)')
    parser.add_argument('-a', '--alternate_subtype_txt', required=True, type=str, default='alternate_subtypes.txt', help='Output TXT file path (will have just the reference IDs)')
    parser.add_argument('--salmon', action='store_true', help='Salmon was used to generate the input (rather than kallisto)')
    parser.add_argument('-o', '--output', type=str, default='best_reference.tsv', help='Output TSV file path')
    parser.add_argument('-t', '--tpm_threshold', type=float, default=250000, help='Threshold for the TPM of reads (quasi-)mapping to a reference for finding other potential subtypes')
    return parser.parse_args()

def process_data(input_path, salmon, sample_name, top_percent, best_ref_txt, alternate_subtype_txt, output_path, fraction_reads):
    df = pd.read_csv(input_path, sep='\t')
    
    # Update headers based on tool used
    if salmon:
        new_headers = ['target_id', 'length', 'eff_length', 'tpm', 'est_counts']
        df.columns = new_headers

    # Choose the best target_id based on the highest tpm
    best_row = df.loc[df['tpm'].idxmax()]
    best_tpm = best_row['tpm']
    best_target_id = best_row['target_id']
    
    # Write the best target id
    with open(best_ref_txt, 'w') as f:
        f.write(f'{best_target_id}\n')
    
    # Extract genotype and subtype
    best_genotype, best_subtype = extract_genotype_subtype(best_target_id)

    # Find the n-best target_ids within top_percent of the highest tpm
    threshold_tpm = best_tpm - (best_tpm * (top_percent / 100))
    close_hits_df = df[df['tpm'] >= threshold_tpm]
    close_hits_list = [x for x in close_hits_df['target_id'].tolist() if x != best_target_id]
    close_hits = ';'.join(close_hits_list)

    # Check for conflicting subtypes where read fraction >= fraction_reads
    other_subtypes = []
    for _, row in df.iterrows():
        if row['target_id'] != best_target_id:
            # read_fraction = row['est_counts'] / total_reads
            read_fraction = row['tpm']
            if read_fraction >= fraction_reads:
                genotype, subtype = extract_genotype_subtype(row['target_id'])
                if genotype+subtype != best_genotype+best_subtype:
                    other_subtypes.append(row['target_id'])

    # also consider the close_hits from earlier
    for ch in close_hits_list:
        genotype, subtype = extract_genotype_subtype(ch)
        if genotype+subtype != best_genotype+best_subtype:
            other_subtypes.append(ch)
    
    # make sure there are no dupes
    other_subtypes = list(set(other_subtypes))
    
    if len(other_subtypes) > 0:
        print(f"Warning: the main subtype was detected as {best_genotype+best_subtype}, but there are potentially other subtypes present: {', '.join(other_subtypes)}")
    
    with open(alternate_subtype_txt, 'w') as f:
        # also write in the main subtype to help with read mapping 
        f.write(f'{best_target_id}\n')
        for st in other_subtypes:
            f.write(f'{st}\n')

    # Write out the TSV file
    with open(output_path, 'w') as f:
        f.write(f'sample_id\tgenotype\tsubtype\tbest_ref\tclose_hits\tother_potential_subtypes\n')
        f.write(f'{sample_name}\t{best_genotype}\t{best_subtype}\t{best_target_id}\t{close_hits}\t{";".join(other_subtypes)}\n')

def extract_genotype_subtype(target_id):
    parts = target_id.split('_')
    genotype = ''.join([char for char in parts[0] if char.isdigit()])
    subtype = parts[0][len(genotype):]
    return genotype, subtype

if __name__ == "__main__":
    args = parse_args()
    process_data(args.input, args.salmon, args.sample_name, args.top, args.best_ref_txt, args.alternate_subtype_txt, args.output, args.tpm_threshold)