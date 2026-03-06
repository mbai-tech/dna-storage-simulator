"""
DNA Storage Cost Model — Heckel et al., arXiv:1705.04732v1

Theoretical rates and cost formulae derived from the paper's main results.

  Rs_theory  = (1 - e^{-c}) * (1 - 1/beta)          [Theorem 1, eq. 3]
  Rr_theory  = Rs_theory / c                          [Corollary 1: Rs = c * Rr]
  cost       = q / Rs + 1 / Rr                        [Section 4]
             = (q + c) / ((1 - e^{-c}) * (1 - 1/beta))

where:
  c    = N / M  coverage depth (reads per unique strand)
  beta = L / log2(M)  strand-length scaling factor
  q    = synthesis cost / sequencing cost  (~10,000-100,000 in practice)
"""

import math


def Rs_theory(c: float, beta: float) -> float:
    """
    Asymptotic storage rate: bits stored per nucleotide synthesized.

    Rs = (1 - e^{-c}) * (1 - 1/beta)

    Returns 0 when beta <= 1 (no positive rate achievable).
    """
    if beta <= 1:
        return 0.0
    return (1 - math.exp(-c)) * (1 - 1 / beta)


def Rr_theory(c: float, beta: float) -> float:
    """
    Asymptotic recovery rate: bits recovered per nucleotide sequenced.

    Rr = Rs / c  (from the relation Rs = c * Rr, Corollary 1)

    Returns 0 when c <= 0.
    """
    if c <= 0:
        return 0.0
    return Rs_theory(c, beta) / c


def cost_theory(c: float, q: float, beta: float) -> float:
    """
    Cost per stored bit = q / Rs + 1 / Rr.

    Equivalent closed form: (q + c) / ((1 - e^{-c}) * (1 - 1/beta))

    Args:
        c    : coverage depth
        q    : synthesis-to-sequencing cost ratio
        beta : strand-length scaling factor

    Returns math.inf when Rs = 0 (beta <= 1).
    """
    rs = Rs_theory(c, beta)
    if rs == 0:
        return math.inf
    rr = Rr_theory(c, beta)
    return q / rs + 1 / rr
