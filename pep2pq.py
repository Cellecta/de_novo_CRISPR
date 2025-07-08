import pandas as pd
from Bio import SeqIO

records = []
with open("Homo_sapiens.GRCh38.pep.all.fa") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # Example header: ENST00000373020.9 pep:known chromosome:GRCh38:13:32316461:32357807:1 gene:ENSG00000120738.17

        transcript_id = record.description.split(' ')[4].split(':')[1]
        
        transcript_base = transcript_id.split('.')[0]
        version = transcript_id.split('.')[1]
        seq = str(record.seq)
        molecule = 'protein'
        desc = 'None'
        id = record.id.split('.')[0]
        aa_len = len(seq)
        transcript_len = aa_len*3
        records.append({"Target Transcript": transcript_id, "Target Total Length": transcript_len, "Transcript Base": transcript_base, "version": version, "seq": seq, "molecule": molecule, "desc": desc, "id": id, "AA len": aa_len})

df = pd.DataFrame(records)
df.to_parquet("aa_seqs.pq", index=False)

