# crispr_pipeline/__main__.py
import argparse
from pathlib import Path
import time
from runner import run_pipeline

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('targets', help='Path to FASTA file with target sequences')
    parser.add_argument('reference', help='Path to FASTA file with reference genome sequences')
    parser.add_argument('gtf', help='Path to GTF file for genome annotation')
    parser.add_argument('save', help='Path to save results (TSV)')
    parser.add_argument('--exclude', default='^$', help='Regex for reference seq IDs to exclude')
    parser.add_argument('--nontargeting', action='store_true', help='Filter out targeting guides')
    parser.add_argument('--number', type=int, default=10, help='Number of guides to pick per target')
    
    args = parser.parse_args()

    start = time.time()
    save_path = Path(args.save)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    gtf_path=args.gtf

    with save_path.open('w') as out_file:
        for row in run_pipeline(
            targets=args.targets,
            reference=args.reference,
            exclude=args.exclude,
            nontargeting=args.nontargeting,
            number=args.number,
            temp_path=save_path.parent,
            gtf_path = gtf_path

        ):
            out_file.write('\t'.join(map(str, row)) + '\n')

    print(f"# Done in {time.time() - start:.2f} seconds")

if __name__ == '__main__':
    main()
