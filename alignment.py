import os
import re
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scoring import calcHitScore, calcCfdScore, calcMitGuideScore
from pathlib import Path
from Bio.Restriction import RestrictionBatch

def check_restriction_sites(dna_seq: str, enzymes=None):
    """
    Check if the input DNA sequence is cut by any restriction enzyme.

    Parameters:
        dna_seq (str): DNA sequence (e.g., 20-mer guide).
        enzymes (list or None): List of enzyme names (str), or None for common enzymes.

    Returns:
        list of enzyme names that cut the sequence.
    """
    seq_obj = Seq(dna_seq)

    if enzymes is None:
        # Define commonly used enzymes
        enzymes = [
            "BsaI", "BsmBI", "EcoRI", "BamHI", "HindIII", "XhoI", "NheI", "NotI", "ApaI", "SmaI", "SpeI", "KpnI"
        ]

    batch = RestrictionBatch(enzymes)
    result = batch.search(seq_obj)

    # Return enzyme names that cut the sequence
    return [enzyme.__name__ for enzyme, sites in result.items() if sites]

def get_pams(seq: str):
    """
    Scan both strands of the input DNA sequence to find PAM sites (NGG or CCN).

    Yields dictionaries with:
        - position: 0-based start of 30-mer window
        - context_30nt: 30 nt context (25 upstream + PAM + 3 downstream)
        - guide_seq: 20 nt sgRNA
        - pam: PAM sequence (NGG)
        - strand: "sense" or "antisense"
        - cut_site: cut site position (0-based index in sequence)
        - cut_site_pct: cut_site / sequence_length * 100 (%)
    """
    seq_len = len(seq)

    # Sense strand: NGG
    for match in re.finditer(r'(?=([ACGT]{25}GG[ACGT]{3}))', seq, re.IGNORECASE):
        context_30nt = match.group(1)
        guide_seq = context_30nt[4:24]
        pam = context_30nt[24:27]
        guide_pam = guide_seq + pam
        position = match.start()
        cut_site = position + 17
        yield {
            "position": position,
            "context_30nt": context_30nt,
            "guide_seq": guide_seq,
            "pam": pam,
            "guide_pam": guide_pam,
            "strand": "sense",
            "cut_site": cut_site,
            "cut_site_pct": round(cut_site / seq_len * 100, 2),
            "restriction_enzymes": check_restriction_sites(guide_seq)

        }

    # Antisense strand: CCN
    for match in re.finditer(r'(?=([ACGT]{3}CC[ACGT]{25}))', seq, re.IGNORECASE):
        context = match.group(1)
        rc_30mer = str(Seq(context).reverse_complement())
        guide_seq = rc_30mer[4:24]
        pam = rc_30mer[24:27]
        guide_pam = guide_seq + pam
        position = match.start()
        cut_site = position + 30 - 17
        yield {
            "position": position,
            "context_30nt": rc_30mer,
            "guide_seq": guide_seq,
            "pam": pam,
            "guide_pam": guide_pam,
            "strand": "antisense",
            "cut_site": cut_site,
            "cut_site_pct": round(cut_site / seq_len * 100, 2),
            "restriction_enzymes": check_restriction_sites(guide_seq)
        }


def run_bwa(refs_path: str, seqs: list[str], temp_path: str):
    """
    Runs BWA on a set of sequences and generates .sam output.
    """
    required_ext = ['amb', 'ann', 'bwt', 'pac', 'sa']
    if not all(Path(f"{refs_path}.{ext}").exists() for ext in required_ext):
        subprocess.run(f'bwa index "{refs_path}"', shell=True, check=True)

    temp_fna = os.path.join(temp_path, "temp.fna")
    with open(temp_fna, 'w') as f:
        f.write('\n'.join(f">{s}\n{s}" for s in seqs))

    temp_sai = os.path.join(temp_path, "temp.sai")
    subprocess.run(
        f'bwa aln -k 4 -l 20 -n 4 -o 0 -N -t 4 "{refs_path}" "{temp_fna}" > "{temp_sai}"',
        shell=True, check=True
    )

    temp_sam = os.path.join(temp_path, "temp.sam")
    subprocess.run(
        f'bwa samse -n 60000 "{refs_path}" "{temp_sai}" "{temp_fna}" > "{temp_sam}"',
        shell=True, check=True
    )
    return temp_sam



