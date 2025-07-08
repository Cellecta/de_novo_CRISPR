# crispr_pipeline/__main__.py
import argparse
from pathlib import Path
import time
from runner import run_pipeline
from offtarget_annotator import build_feature_trees_cached, build_feature_trees
import pandas as pd
from Bio import SeqIO
#from scoring import rs3_target_score

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('targets', help='Path to FASTA file with target sequences')
    parser.add_argument('reference', help='Path to FASTA file with reference genome sequences')
    parser.add_argument('gtf', help='Path to GTF file for genome annotation')
    parser.add_argument('save', help='Path to save results (TSV)')
    parser.add_argument('--exclude', default='^$', help='Regex for reference seq IDs to exclude')
    parser.add_argument('--nontargeting', action='store_true', help='Filter out targeting guides')
    parser.add_argument('--number', type=int, default=10, help='Number of guides to pick per target')
    parser.add_argument('--transcript_anno', required=True, help='Path to gene_mane_summary.tsv file')
    
    args = parser.parse_args()

    start = time.time()
    save_path = Path(args.save)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    gtf_path=args.gtf


    print("[INFO] Loading gene annotations...")
    gtf_trees = build_feature_trees_cached(gtf_path)


# Get header and data separately
    rows = list(run_pipeline(
        targets=args.targets,
        reference = args.reference,
        exclude=args.exclude,
        nontargeting=args.nontargeting,
        number=args.number,
        temp_path=save_path.parent,
        gtf_trees = gtf_trees,
        transcript_anno=args.transcript_anno
    ))

    header, *data = rows
    df = pd.DataFrame(data, columns=header)

    # Save output
    df.to_csv(save_path, sep='\t', index=False)

    print(f"✅ Finished. Saved {len(df)} rows to {save_path} in {time.time() - start:.2f} seconds.")



if __name__ == '__main__':
    main()
