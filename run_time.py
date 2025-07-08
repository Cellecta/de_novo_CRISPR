from Bio import SeqIO
import time

start = time.time()

fasta_file = "/run/media/Dongfang/DataUSB/Cellecta_Project/denovo_crispr/db/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"  # or hg19.fa
records = list(SeqIO.parse(fasta_file, "fasta"))

end = time.time()
print(f"Loaded {len(records)} sequences in {end - start:.2f} seconds")

