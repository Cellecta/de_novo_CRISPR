import os
import pickle
from pathlib import Path
from typing import Tuple

# Global cache for scores
_mm_scores, _pam_scores = None, None

def get_data_path(filename: str) -> Path:
    """Return the full path to a data file in the current package directory."""
    return Path(__file__).parent / filename

def load_pickle_data(filepath: Path):
    """Load a pickle file and raise informative error if missing."""
    if not filepath.exists():
        raise FileNotFoundError(f"Required data file not found: {filepath}")
    with open(filepath, 'rb') as f:
        return pickle.load(f)

def get_mm_pam_scores() -> Tuple[dict, dict]:
    """Load mismatch and PAM score dictionaries (cached)."""
    global _mm_scores, _pam_scores
    if _mm_scores is None or _pam_scores is None:
        base_dir = Path(__file__).parent
        _mm_scores = load_pickle_data(base_dir / 'mismatch_score.pkl')
        _pam_scores = load_pickle_data(base_dir / 'pam_scores.pkl')
    return _mm_scores, _pam_scores

def revcom(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def ensure_dir(path: Path):
    """Create directory if it doesn’t exist."""
    path.mkdir(parents=True, exist_ok=True)

def read_fasta_ids(fasta_path: str) -> list:
    """Read only FASTA sequence IDs."""
    from Bio import SeqIO
    return [record.id for record in SeqIO.parse(fasta_path, 'fasta')]
