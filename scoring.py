import os
import re
import pickle
from utils import get_mm_pam_scores, revcom, get_data_path, load_pickle_data

HIT_SCORE_WEIGHTS = [
    0.000, 0.000, 0.014, 0.000, 0.000, 0.395, 0.317, 0.000, 0.389, 0.079,
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583
]

def calcHitScore(guide: str, target: str) -> float:
    """MIT off-target score"""
    if len(guide) == len(target) == 23:
        guide = guide[:20]
        target = target[:20]

    mismatches = []
    score = 1.0
    last_mismatch = None

    for i in range(20):
        if guide[i] != target[i]:
            mismatches.append(i)
            if last_mismatch is not None:
                pass  # distance not used directly here
            score *= 1 - HIT_SCORE_WEIGHTS[i]
            last_mismatch = i

    mm_count = len(mismatches)

    if mm_count < 2:
        dist_penalty = 1.0
    else:
        avg_dist = sum(b - a for a, b in zip(mismatches, mismatches[1:])) / len(mismatches[1:])
        dist_penalty = 1.0 / (((19 - avg_dist) / 19.0) * 4 + 1)

    count_penalty = 1.0 if mm_count == 0 else 1.0 / (mm_count ** 2)

    return score * dist_penalty * count_penalty * 100

def calcMitGuideScore(score_sum: float) -> int:
    score = 100 / (100 + score_sum)
    return int(round(score * 100))

def revcom(seq: str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    return ''.join(complement[base] for base in reversed(seq))

# --- CFD Score Loading ---

_mm_scores, _pam_scores = None, None

def get_mm_pam_scores():
    global _mm_scores, _pam_scores
    if _mm_scores is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        _mm_scores = load_pickle_data(get_data_path('mismatch_score.pkl'))
        _pam_scores = load_pickle_data(get_data_path('pam_scores.pkl'))
    return _mm_scores, _pam_scores

def calc_cfd(wt: str, sg: str, pam: str) -> float:
    mm_scores, pam_scores = get_mm_pam_scores()
    score = 1.0
    sg = sg.replace('T', 'U')
    wt = wt.replace('T', 'U')

    for i, (wt_base, sg_base) in enumerate(zip(wt, sg)):
        if wt_base != sg_base:
            key = f"r{wt_base}:d{revcom(sg_base)},{i+1}"
            score *= mm_scores.get(key, 1.0)

    score *= pam_scores.get(pam, 0.0)
    return score

def calcCfdScore(guide: str, off_target: str) -> float:
    guide = guide.upper()
    off_target = off_target.upper()
    if not re.search('[^ATCG]', guide + off_target):
        pam = off_target[-2:]
        sg = off_target[:20]
        return calc_cfd(guide, sg, pam)
    return 0.0


# --- RuleSet3 + Target ---









