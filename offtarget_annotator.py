
import argparse
import json
from Bio import SeqIO
from gtfparse import read_gtf
from intervaltree import IntervalTree
from collections import defaultdict
from scoring import calcCfdScore, calcMitGuideScore  # Make sure scoring.py is available
import pickle
import os
from Bio.Seq import Seq

def parse_xa_tag(xa_string):
    xa_data = xa_string.replace('XA:Z:', '').strip().split(';')
    for entry in xa_data:
        if not entry:
            continue
        chrom, strand_pos, cigar, mm = entry.split(',')
        strand = strand_pos[0]
        pos = int(strand_pos[1:])
        yield chrom, strand, pos, cigar, int(mm)


        

def build_feature_trees(gtf_file):
    gtf = read_gtf(gtf_file)
    trees = {}
    exon_map = defaultdict(list)

    for _, row in gtf.iterrows():
        feature = row['feature']
        if feature not in {"exon", "CDS", "UTR", "gene"}:
            continue

        chrom = row['seqname']
        start, end = int(row['start']) - 1, int(row['end'])
        meta = {
            'feature': feature,
            'gene_id': row.get('gene_id', ''),
            'transcript_id': row.get('transcript_id', ''),
            'gene_type': row.get('gene_type') or row.get('gene_biotype', ''),
            'gene_name': row.get('gene_name', ''),
            'strand': row.get('strand', '.')
        }

        trees.setdefault(chrom, {}).setdefault(feature, IntervalTree())[start:end] = meta

        if feature == "exon" and meta['transcript_id']:
            exon_map[(chrom, meta['transcript_id'])].append((start, end))

    for (chrom, tid), exons in exon_map.items():
        exons.sort()
        for (end1, start2) in zip([e[1] for e in exons[:-1]], [e[0] for e in exons[1:]]):
            if start2 > end1:
                meta = {
                    'feature': 'intron',
                    'gene_id': '',
                    'transcript_id': tid,
                    'gene_type': '',
                    'gene_name': '',
                    'strand': '.'
                }
                trees.setdefault(chrom, {}).setdefault('intron', IntervalTree())[end1:start2] = meta

    return trees

def build_feature_trees_cached(gtf_file, cache_file=None):
    if cache_file is None:
        cache_file = gtf_file + ".tree.pkl"

    if os.path.exists(cache_file):
        with open(cache_file, 'rb') as f:
            return pickle.load(f)

    # Build from scratch
    trees = build_feature_trees(gtf_file)

    with open(cache_file, 'wb') as f:
        pickle.dump(trees, f)

    return trees

def annotate_hit(chrom, pos, strand, trees, window=23):
    annotations = []
    cds_genes = set()
    start = pos - 1
    end = start + window
    if chrom not in trees:
        return ["Intergenic"], []

    for feature in ["CDS", "UTR", "exon", "intron"]:
        if feature not in trees[chrom]:
            continue
        for hit in trees[chrom][feature][start:end]:
            gid = hit.data.get('gene_id', '')
            tid = hit.data.get('transcript_id', '')
            gtype = hit.data.get('gene_type', '')
            gname = hit.data.get('gene_name', '')
            label = f"{feature}:{gid}/{tid} [{gtype}, {gname}]" if gname else f"{feature}:{gid}/{tid} [{gtype}]"
            annotations.append(label)
            if feature == "CDS" and gid:
                cds_genes.add((gid, gname))

    cds_gene_list = [{"gene_id": gid, "gene_name": gname} for gid, gname in sorted(cds_genes)]
    return annotations or ["Intergenic"], cds_gene_list

def extract_sequence(genome_dict, chrom, strand, pos, length=23):
    start = pos - 1
    seq = genome_dict[chrom].seq[start:start+length]
    return str(seq.reverse_complement() if strand == '-' else seq)

def summarize_results_by_guide(results, summary_output=None):
    summary_dict = {}

    for guide_seq, alignments in results.items():
        total = len(alignments)
        match_1_20 = 0
        match_3_20 = 0
        match_1_20_cds = 0
        match_3_20_cds = 0
        cfd_all = []
        cfd_cds = []

        for aln in alignments:
            cfd = aln.get("cfd_score", 0.0)
            cds_hit = aln.get("cds_gene_count", 0) > 0
            cfd_all.append(cfd)
            if cds_hit:
                cfd_cds.append(cfd)

            if aln.get("perfect_match_1_20"):
                match_1_20 += 1
                if cds_hit:
                    match_1_20_cds += 1

            if aln.get("perfect_match_3_20"):
                match_3_20 += 1
                if cds_hit:
                    match_3_20_cds += 1

        avg_cfd = round(calcMitGuideScore(sum(cfd_all)) if cfd_all else 0.0, 4)
        avg_cfd_cds = round(calcMitGuideScore(sum(cfd_cds)) if cfd_cds else 0.0, 4)

        summary_dict[guide_seq] = {
            "total_alignments": total,
            "perfect_match_1_20": match_1_20,
            "perfect_match_3_20": match_3_20,
            "perfect_match_1_20_CDS": match_1_20_cds,
            "perfect_match_3_20_CDS": match_3_20_cds,
            "cfd_all": avg_cfd,
            "cfd_CDS": avg_cfd_cds
        }

    # Optional: save to CSV or JSON
    if summary_output:
        if summary_output.endswith(".json"):
            with open(summary_output, "w") as f:
                json.dump(summary_dict, f, indent=2)
        elif summary_output.endswith(".csv"):
            import csv
            # flatten dict for CSV output; each row gets the guide as a column
            fieldnames = ["guide", "total_alignments", "perfect_match_1_20",
                          "perfect_match_3_20", "perfect_match_1_20_CDS",
                          "perfect_match_3_20_CDS", "cfd_all", "cfd_CDS"]
            with open(summary_output, "w", newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for guide, data in summary_dict.items():
                    row = {"guide": guide}
                    row.update(data)
                    writer.writerow(row)

    return summary_dict






def main(sam_file, genome_fasta, gtf_path):
    print("[INFO] Loading reference genome...")
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    print("[INFO] Loading gene annotations...")
    gtf_trees = build_feature_trees_cached(gtf_path)

    results = defaultdict(list)

    print(f"[INFO] Reading SAM file: {sam_file}")
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            guide_seq = fields[0].upper()
            flag = int(fields[1])
            chrom = fields[2]
            pos = int(fields[3])
            strand = '-' if (flag & 16) else '+'
            annots, cds_gene_list = annotate_hit(chrom, pos, strand, gtf_trees)

            results[guide_seq].append({
                    "chrom": chrom,
                    "strand": strand,
                    "position": pos,
                    "mismatches": 0,
                    "sequence": guide_seq,
                    "annotations": annots,
                    "cfd_score": 1.0,
                    "perfect_match_1_20": True,
                    "perfect_match_3_20": True,
                    "cds_genes": cds_gene_list,
                    "cds_gene_count": len(cds_gene_list)
                })



            xa_field = next((f for f in fields if f.startswith("XA:Z:")), None)
            if not xa_field:
                continue

            for chrom, strand, pos, cigar, mm in parse_xa_tag(xa_field):
                if chrom not in genome:
                    continue

                seq =extract_sequence(genome, chrom, strand, pos)
                annots, cds_gene_list = annotate_hit(chrom, pos, strand, gtf_trees)
                cfd = calcCfdScore(guide_seq, seq)
                seq = seq.upper()
                pam_seq = seq[-3:]
                target_seq = seq[:20]


                has_ngg_pam = pam_seq[1:] == 'GG' and pam_seq[0] in 'ACGT'
                perfect_1_20 = guide_seq[:20] == target_seq if has_ngg_pam else False
                perfect_3_20 = guide_seq[2:20] == target_seq[2:20] if has_ngg_pam else False

                if has_ngg_pam:
                    annots.append('NGG_PAM')
                if perfect_1_20:
                    annots.append("Match1-20")
                if perfect_3_20:
                    annots.append("Match3-20")

                results[guide_seq].append({
                    "chrom": chrom,
                    "strand": strand,
                    "position": pos,
                    "mismatches": mm,
                    "sequence": seq,
                    "annotations": annots,
                    "cfd_score": round(cfd, 4) if cfd else 0.0,
                    "perfect_match_1_20": perfect_1_20,
                    "perfect_match_3_20": perfect_3_20,
                    "cds_genes": cds_gene_list,
                    "cds_gene_count": len(cds_gene_list)
                })

    return summarize_results_by_guide(results)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", default="summary.json")
    args = parser.parse_args()

    result = main(args.sam, args.fasta, args.gtf)
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)
