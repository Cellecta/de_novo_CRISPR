from Bio import SeqIO
import numpy as np
from design import design_test
from scoring import calcMitGuideScore
from alignment import get_pams, run_bwa
from rs3.seq import predict_seq
import pandas as pd
from offtarget_annotator import main as run_annotation
import json


def load_transcript_annotations(transcript_anno):
    df = pd.read_csv(transcript_anno, sep="\t")
    df['Transcript_id'] = df['Transcript_ID'].str.split(',')
    df = df.explode('Transcript_id')
    df['Transcript_id'] = df['Transcript_id'].str.replace(r'\.\d+$', '', regex=True)
    return df.set_index('Transcript_id').to_dict(orient='index')



def run_pipeline(targets, reference, exclude, nontargeting, number, temp_path, gtf_trees, transcript_anno):

    yield [
    'Index', 'Transcript_ID', 'Gene_ID', 'Gene_Symbol', 'Isoform_ID',
    'guide', 'PAM', 'Context_30nt', 'Orientation', 'Chromosome', 'Strand',
    'Cutting Position', 'Target Cut %', 'Restriction_Enzyme', 'RS3_score',
    'perfect_match_1_20_all', 'perfect_match_1_20_cds',
    'perfect_match_3_20_all', 'perfect_match_3_20_cds',
    'CFD_score_all', 'CFD_score_cds', 'Isoform_fraction', 'gnomAD'
    ]

    transcript_anno = load_transcript_annotations(transcript_anno)
    



    for n, record in enumerate(SeqIO.parse(targets, 'fasta')):
        print(f'# Target={record.id} length={len(record.seq)} bp')
        hits = list(get_pams(str(record.seq)))
        cds_id_list = [str(record.id).split('.')[0]]*len(hits)
        full_id = str(record.id)
        base_id = record.id.split('.')[0]
        # Fetch annotation info
        gene_info = transcript_anno.get(base_id, {})
        gene_id = gene_info.get('Gene_ID', '')
        gene_symbol = gene_info.get('Gene_Symbol', '')
        refseq_ids = gene_info.get('Transcript_ID', '')
        refseq_ids_list = refseq_ids.split(',') if isinstance(refseq_ids, str) else [refseq_ids]

        if not hits:
            continue
        

        rs3_seqs = [h["context_30nt"] for h in hits]
        rs3_scores = predict_seq(rs3_seqs, sequence_tracr='Chen2013')

        guide_pams = [h["guide_pam"] for h in hits]

        # modified parse_bwa that also returns the SAM path
        run_bwa( reference, guide_pams, temp_path)

        # run annotation script from sam_path
        offtarget_dict = run_annotation(temp_path/'temp.sam', reference, gtf_trees, cds_id_list, refseq_ids_list)

        results = []
        for i, hit in enumerate(hits):
            guide = hit["guide_seq"]
            guide_pam = hit["guide_pam"]
            context_30nt = hit['context_30nt']
            rs3_score = rs3_scores[i]
            if rs3_score is None:
                continue
            

            # Get all per-site entries for this guide sequence
            matched_keys = [key for key in offtarget_dict if key.startswith(f"{guide_pam}_")]

            if not matched_keys:
                # No off-target info, still yield minimal info
                results.append((
                    guide, hit["pam"], hit['strand'], None, None, None, hit["cut_site_pct"], hit['restriction_enzymes'],
                    rs3_score, 0, 0, 0, 0, 0.0, 0.0, 0.0, ''
                ))
                continue

            for key in matched_keys:

                stats = offtarget_dict.get(key, {})

                p1_all = stats.get('perfect_match_1_20', 0)
                p1_cds = stats.get('perfect_match_1_20_CDS', 0)
                p3_all = stats.get('perfect_match_3_20', 0)
                p3_cds = stats.get('perfect_match_3_20_CDS', 0)
                cfd_all = stats.get('cfd_all', 0.0)
                cfd_cds = stats.get('cfd_CDS', 0.0)
                chr = stats.get('chr', 0.0)
                cut_pos = stats.get('cut_pos', 0.0)
                strand = stats.get('strand', 0.0)
                Isoform_fra = stats.get('Isoform_fra', 0.0)
                gnomAD = stats.get('gnomAD', 0.0)

                results.append((
                    guide, hit["pam"],context_30nt, hit['strand'], chr, strand, cut_pos, hit["cut_site_pct"], hit['restriction_enzymes'],
                    rs3_score, p1_all, p1_cds, p3_all, p3_cds,
                    round(cfd_all, 4), round(cfd_cds, 4), Isoform_fra, gnomAD
                ))

        results.sort(key=lambda x: x[9], reverse=True)
        for row in results[:number]:
            yield (n + 1, full_id, gene_id, gene_symbol, refseq_ids) + row
