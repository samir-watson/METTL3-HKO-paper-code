# environment= transcript2genome
# Before using a specific Ensembl release for the first time, you need to download and index its data. You can do this using the pyensembl install command in your terminal
# pyensembl install --release 99 --species human

import os
import pandas as pd
from pyensembl import EnsemblRelease
import pybedtools

# Specify the Ensembl release
data = EnsemblRelease(99)

# Set the base directory
base_dir = '/home/samirwatson/faststorage/METTL3/DRS_data/xpore/diff_mod/all_KO/'

# Function to convert Ensembl chromosome format to GENCODE format
def convert_to_gencode(chromosome):
    if chromosome == "MT":
        return "chrM"
    else:
        return "chr" + chromosome

def transcript_to_genome(id, position):
    if id.strand == '+':
        return id.start + position - 1
    else:
        return id.end - position + 1

# Function to convert transcript coordinates to genomic coordinates

# Find all CSV files in the base directory and its subdirectories that match the pattern
csv_files = [os.path.join(root, f) for root, dirs, files in os.walk(base_dir) for f in files if f.startswith('all_KO_filtered_majority_direction_kmer_diffmod') and f.endswith('.7.table')]

malformed_entries = 0

for csv_file in csv_files:
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file)
    # Add columns for chromosome and strand
    df['chromosome'] = df['id'].apply(lambda x: convert_to_gencode(data.transcript_by_id(x.split('.')[0]).contig) if data.transcript_by_id(x.split('.')[0]) else None)
    df['strand'] = df['id'].apply(lambda x: data.transcript_by_id(x.split('.')[0]).strand if data.transcript_by_id(x.split('.')[0]) else None)
    # Calculate start and stop positions before calling transcript_to_genome
    df['start_pos'] = df.apply(lambda row: transcript_to_genome(data.transcript_by_id(row['id'].split('.')[0]), row['position'] - 2) if data.transcript_by_id(row['id'].split('.')[0]) else None, axis=1)
    df['end_pos'] = df.apply(lambda row: transcript_to_genome(data.transcript_by_id(row['id'].split('.')[0]), row['position'] + 2) if data.transcript_by_id(row['id'].split('.')[0]) else None, axis=1)
    # Swap start and end positions if transcript is on the negative strand
    df.loc[df['strand'] == '-', ['start_pos', 'end_pos']] = df.loc[df['strand'] == '-', ['end_pos', 'start_pos']].values
    # Check for malformed entries
    malformed = df[df['start_pos'] > df['end_pos']]
    malformed_entries += len(malformed)
    # Merge transcript_id with kmer
    df['id'] = df['id'] + '_' + df['kmer']
    # Create a BED file from the DataFrame
    bed_df = df[['chromosome', 'start_pos', 'end_pos', 'id', 'pval_WT_vs_xKO', 'strand']]
    bed_df.to_csv(os.path.join(os.path.dirname(csv_file), f'{os.path.basename(csv_file)}_transcript_positions.bed'), sep='\t', index=False, header=False)
    transcript_bed = pybedtools.BedTool(os.path.join(os.path.dirname(csv_file), f'{os.path.basename(csv_file)}_transcript_positions.bed'))
    # Sort the transcript BED file
    sorted_transcript_bed = transcript_bed.sort()

print(f"Number of malformed entries: {malformed_entries}")
