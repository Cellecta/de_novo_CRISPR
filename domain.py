import pandas as pd
from rs3.targetdata import write_transcript_data


# Load your aa_seqs.pq file
design_df = pd.read_parquet("aa_seqs.pq")

# Make sure 'transcript_id' and 'aa_sequence' columns exist
# Add required 'Target Total Length' column
design_df['Target Total Length'] = design_df['aa_sequence'].str.len()



# Run write_transcript_data
write_transcript_data(
    design_df=design_df,
    save_file="domains.pq",
    transcript_id_col="transcript_id",
    transcript_len_col="Target Total Length"
)
