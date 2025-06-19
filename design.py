import re
from Bio.Restriction import BsaI
from Bio.Seq import Seq

# Predefined sgRNA design templates
DESIGNS = [((BsaI,), 'TGTGGAAAGGACGAAACACCG{}GTTTAAGAGCTATGCTGGAAGGAGA')]

def design_test(design_index: int):
    """
    Returns a function that checks whether an sgRNA is cloneable.
    It must not:
      - contain GGGGG at the start or TTTT anywhere
      - contain a BsaI restriction site in the full oligo
    """
    enzymes, template = DESIGNS[design_index]
    pattern = re.compile(r'(^[Gg]{5})|([Tt]{4})')

    def _is_cloneable(seq: str) -> bool:
        full_seq = Seq(template.format(seq))
        has_restriction_site = any(enzyme.search(full_seq) for enzyme in enzymes)
        has_bad_motif = bool(pattern.search(seq))
        return not (has_bad_motif or has_restriction_site)

    return _is_cloneable
