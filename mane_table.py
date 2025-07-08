import pandas as pd

# Load MANE summary file
mane = pd.read_csv("MANE.GRCh38.v1.4.summary.txt", sep="\t")
mane.columns = [col.lstrip('#') for col in mane.columns]  # clean column headers

# Optional: remove Ensembl version numbers
mane['Ensembl_Gene'] = mane['Ensembl_Gene'].str.replace(r'\.\d+$', '', regex=True)
mane['Ensembl_nuc'] = mane['Ensembl_nuc'].str.replace(r'\.\d+$', '', regex=True)
mane['RefSeq_nuc'] = mane['RefSeq_nuc'].str.replace(r'\.\d+$', '', regex=True)

# Group and aggregate
grouped = (
    mane.groupby(['Ensembl_Gene', 'symbol'])
    .agg({
        'NCBI_GeneID': 'first',
        'HGNC_ID': 'first',
        'name': 'first',
        'RefSeq_nuc': lambda x: ','.join(sorted(set(x))),
        'Ensembl_nuc': lambda x: ','.join(sorted(set(x))),
        'MANE_status': lambda x: ','.join(sorted(set(x))),
        'GRCh38_chr': 'first',
        'chr_start': 'first',
        'chr_end': 'first',
        'chr_strand': 'first'
    })
    .reset_index()
)

# Save to file
grouped.to_csv("gene_mane_summary.tsv", sep="\t", index=False)

# Preview
print(grouped.head())



