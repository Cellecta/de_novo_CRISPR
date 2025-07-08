import pandas as pd
import re
from collections import defaultdict

def extract_ccds_coding_transcripts_with_position(gtf_file):
    cds_transcripts = set()
    gene_info = defaultdict(lambda: {
        "gene_id": "", 
        "gene_name": "", 
        "transcripts": set(),
        "ccds_ids": set(),
        "chr": "", 
        "start": float("inf"), 
        "end": 0, 
        "strand": "",
        "gene_biotype": ""
    })

    # First pass: collect transcripts with CDS
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#') or '\tCDS\t' not in line:
                continue
            tx_match = re.search(r'transcript_id "([^"]+)"', line)
            if tx_match:
                cds_transcripts.add(tx_match.group(1))

    # Second pass: collect transcript info + gene positions and biotype
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue

            feature_type = fields[2]
            chrom = fields[0].removeprefix("chr") if fields[0].startswith("chr") else fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attr = fields[8]

            gene_id_match = re.search(r'gene_id "([^"]+)"', attr)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attr)
            tx_id_match = re.search(r'transcript_id "([^"]+)"', attr)
            biotype_match = re.search(r'transcript_biotype "([^"]+)"', attr)
            gene_biotype_match = re.search(r'gene_biotype "([^"]+)"', attr)
            ccds_match = re.search(r'ccds_id "([^"]+)"', attr)

            if not gene_id_match:
                continue
            gene_id = gene_id_match.group(1)

            if feature_type == "transcript":
                if not (tx_id_match and biotype_match and ccds_match):
                    continue
                tx_id = tx_id_match.group(1)
                transcript_biotype = biotype_match.group(1)

                if transcript_biotype != "protein_coding" or tx_id not in cds_transcripts:
                    continue

                gene_name = gene_name_match.group(1).strip() if gene_name_match and gene_name_match.group(1).strip() else f"novel_protein_{gene_id}"
                ccds_id = ccds_match.group(1)

                gene_info[gene_id]["gene_id"] = gene_id
                gene_info[gene_id]["gene_name"] = gene_name
                gene_info[gene_id]["transcripts"].add(tx_id)
                gene_info[gene_id]["ccds_ids"].add(ccds_id)

            elif feature_type == "gene":
                if gene_info[gene_id]["start"] == float("inf"):
                    gene_info[gene_id]["chr"] = chrom
                    gene_info[gene_id]["start"] = start
                    gene_info[gene_id]["end"] = end
                    gene_info[gene_id]["strand"] = strand
                    if gene_biotype_match:
                        gene_info[gene_id]["gene_biotype"] = gene_biotype_match.group(1)

    # Format final table
    records = []
    for info in gene_info.values():
        if not info["transcripts"]:
            continue
        records.append({
            "Gene_ID": info["gene_id"],
            "Gene_Symbol": info["gene_name"],
            "Chromosome": info["chr"],
            "Start": int(info["start"]) if info["start"] != float("inf") else "",
            "End": info["end"],
            "Strand": info["strand"],
            "Gene_Biotype": info["gene_biotype"],
            "Transcript_ID": ",".join(sorted(info["transcripts"])),
            "CCDS_IDs": ",".join(sorted(info["ccds_ids"]))
        })

    return pd.DataFrame(records)

# Usage
if __name__ == "__main__":
    gtf_path = "Homo_sapiens.GRCh38.114.gtf"  # Update this path as needed
    output_path = "ccds_coding_gene_positions.tsv"

    df = extract_ccds_coding_transcripts_with_position(gtf_path)
    df.to_csv(output_path, sep="\t", index=False)
