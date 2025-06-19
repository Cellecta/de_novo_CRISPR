from Bio import SeqIO
import numpy as np
from design import design_test
from scoring import calcMitGuideScore
from alignment import get_pams, run_bwa
from rs3.seq import predict_seq
from offtarget_annotator import main as run_annotation

def run_pipeline(targets, reference, exclude, nontargeting, number, temp_path, gtf_path):

    yield [
        'Index', 'Transcript_ID', 'Target Cut %', 'strand', 'guide', 'PAM', 'Restriction_Enzyme',
        'RS3_score',
        'perfect_match_1_20_all', 'perfect_match_1_20_cds',
        'perfect_match_3_20_all', 'perfect_match_3_20_cds',
        'CFD_score_all', 'CFD_score_cds'
    ]

    for n, record in enumerate(SeqIO.parse(targets, 'fasta')):
        print(f'# Target={record.id} length={len(record.seq)} bp')
        hits = list(get_pams(str(record.seq)))
        if not hits:
            continue

        rs3_seqs = [h["context_30nt"] for h in hits]
        rs3_scores = predict_seq(rs3_seqs, sequence_tracr='Chen2013')

        guide_pams = [h["guide_pam"] for h in hits]

        # modified parse_bwa that also returns the SAM path
        run_bwa( reference, guide_pams, temp_path)

        # run annotation script from sam_path
        offtarget_dict = run_annotation(temp_path/'temp.sam', reference, gtf_path)

        results = []
        for i, hit in enumerate(hits):
            guide = hit["guide_seq"]
            guide_pam = hit["guide_pam"]
            rs3_score = rs3_scores[i]
            if rs3_score is None:
                continue
            
            stats = offtarget_dict.get(guide_pam, {})
            p1_all = stats.get('perfect_match_1_20', 0)
            p1_cds = stats.get('perfect_match_1_20_CDS', 0)
            p3_all = stats.get('perfect_match_3_20', 0)
            p3_cds = stats.get('perfect_match_3_20_CDS', 0)
            cfd_all = stats.get('cfd_all', 0.0)
            cfd_cds = stats.get('cfd_CDS', 0.0)

            results.append((
                hit["cut_site_pct"], hit["strand"], guide, hit["pam"], hit['restriction_enzymes'],
                rs3_score, p1_all, p1_cds, p3_all, p3_cds,
                round(cfd_all, 4), round(cfd_cds, 4)
            ))

        results.sort(key=lambda x: x[5], reverse=True)
        for row in results[:number]:
            yield (n + 1, record.id) + row
