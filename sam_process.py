import json
from Bio import SeqIO
from gtfparse import read_gtf
from intervaltree import IntervalTree
from collections import defaultdict

# --- 1. Parse XA:Z tag ---
def parse_xa_tag(xa_string: str):
    xa_data = xa_string.replace('XA:Z:', '').strip().split(';')
    results = []
    for entry in xa_data:
        if not entry:
            continue
        chrom, strand_pos, cigar, mismatches = entry.split(',')
        strand = strand_pos[0]
        pos = int(strand_pos[1:])
        results.append((chrom, strand, pos, cigar, int(mismatches)))
    return results

# --- 2. Build genomic annotation trees from GTF ---
def build_feature_trees(gtf_file):
    gtf = read_gtf(gtf_file)
    trees = defaultdict(lambda: defaultdict(IntervalTree))
    for _, row in gtf.iterrows():
        if row['feature'] not in {"exon", "CDS", "UTR", "gene"}:
            continue
        chrom = row['seqname']
        start, end = int(row['start']) - 1, int(row['end'])  # Convert to 0-based
        feature = row['feature']
        meta = {
            'feature': feature,
            'gene_id': row.get('gene_id', ''),
            'transcript_id': row.get('transcript_id', ''),
            'strand': row.get('strand', '.')
        }
        trees[chrom][feature][start:end] = meta
    return trees

# --- 3. Annotate one off-target position ---
def annotate_hit(chrom, pos, strand, trees, window=23):
    annotations = []
    start = pos - 1
    end = start + window
    if chrom not in trees:
        return ["Intergenic"]
    for feature_type, tree in trees[chrom].items():
        for hit in tree[start:end]:
            meta = hit.data
            if strand == meta.get('strand'):
                label = f"{meta['feature']}:{meta['gene_id']}/{meta['transcript_id']}"
                annotations.append(label)
    return annotations or ["Intergenic"]

# --- 4. Extract sequence from genome ---
def extract_sequence(genome_dict, chrom, strand, pos, length=23):
    start = pos - 1
    seq = genome_dict[chrom].seq[start:start + length]
    return str(seq.reverse_complement() if strand == '-' else seq)

# --- 5. Main wrapper: from SAM line to JSON ---
def process_sam_entry(sam_line, genome_fasta, gtf_path, output_json):
    fields = sam_line.strip().split('\t')
    guide_seq = fields[0]
    xa_tag = next((f for f in fields if f.startswith("XA:Z:")), None)
    if not xa_tag:
        print(f"No XA tag for {guide_seq}")
        return

    # Load genome FASTA into memory
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    # Build annotation trees
    trees = build_feature_trees(gtf_path)

    off_target_entries = []
    for chrom, strand, pos, cigar, mm in parse_xa_tag(xa_tag):
        if chrom not in genome:
            continue
        sequence = extract_sequence(genome, chrom, strand, pos)
        annotations = annotate_hit(chrom, pos, strand, trees)
        off_target_entries.append({
            "chrom": chrom,
            "strand": strand,
            "position": pos,
            "mismatches": mm,
            "sequence": sequence,
            "annotations": annotations
        })

    # Save as JSON
    with open(output_json, "w") as f:
        json.dump({guide_seq: off_target_entries}, f, indent=2)
    print(f"Saved output to {output_json}")
