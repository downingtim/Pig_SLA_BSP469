import pandas as pd
import glob
import re

# --- Filtering Parameters ---
# Adjust MIN_LENGTH based on your expected amplicon size (e.g., 400 drops the 185bp fragments)
MIN_LENGTH = 400   
MIN_PIDENT = 98.0  
# Minimum bitscore difference required between the 1st and 2nd best hit to confidently assign an allele
TIE_MARGIN = 2.0   

columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
           "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

all_samples_data = []

for blast_file in glob.glob("*_blast.txt"):
    sample_name = blast_file.replace("_blast.txt", "")
    
    # Load data
    df = pd.read_csv(blast_file, sep='\t', names=columns)
    
    # 1. Strip spurious short alignments and low-identity hits
    df = df[(df['length'] >= MIN_LENGTH) & (df['pident'] >= MIN_PIDENT)]
    
    if df.empty:
        continue
        
    # Extract the read count from the vsearch header (e.g., "size=29484")
    df['read_count'] = df['qseqid'].apply(lambda x: int(re.search(r'size=(\d+)', x).group(1)))
    
    # Sort by read ID, then by highest bitscore
    df = df.sort_values(by=['qseqid', 'bitscore'], ascending=[True, False])
    
    unambiguous_hits = []
    
    # 2. Resolve ambiguous multi-mappers
    for qseqid, group in df.groupby('qseqid'):
        if len(group) == 1:
            unambiguous_hits.append(group.iloc[0])
        else:
            top_hit = group.iloc[0]
            second_hit = group.iloc[1]
            
            # Keep the read only if the top hit beats the second hit by the defined margin
            if (top_hit['bitscore'] - second_hit['bitscore']) >= TIE_MARGIN:
                unambiguous_hits.append(top_hit)

    if not unambiguous_hits:
        continue
        
    filtered_df = pd.DataFrame(unambiguous_hits)
    
    # 3. Sum read depths per confirmed allele for this sample
    allele_counts = filtered_df.groupby('sseqid')['read_count'].sum().reset_index()
    allele_counts['sample'] = sample_name
    all_samples_data.append(allele_counts)

# Combine and pivot into a final Sample x Allele matrix
if all_samples_data:
    master_df = pd.concat(all_samples_data)
    count_matrix = master_df.pivot(index='sample', columns='sseqid', values='read_count').fillna(0)
    
    # Ensure integers for read counts
    count_matrix = count_matrix.astype(int)
    
    count_matrix.to_csv("Filtered_SLA_Allele_Matrix.csv")
    print(f"Matrix generated successfully with {len(count_matrix)} samples.")
else:
    print("No hits passed the filtering thresholds.")
