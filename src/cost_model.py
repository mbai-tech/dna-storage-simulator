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
import numpy as np


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


def optimal_c_theory(q: float, beta: float,
                     c_min: float = 0.01, c_max: float = 50.0,
                     n_points: int = 10000) -> float:
    """
    Find the coverage c* that minimises cost_theory(c, q, beta) analytically.

    The cost simplifies to (q + c) / (1 - e^{-c}) up to the constant factor
    (1 - 1/beta), so c* depends only on q (not on beta or M).  There is no
    closed-form solution, so we use a fine numerical grid over [c_min, c_max].

    The optimality condition is  e^c = q + c + 1  (from d/dc cost = 0), which
    has a unique solution for q > 0.  For q=10,000 this gives c* ≈ 9.2,
    matching the paper (Section 4).

    Args:
        q        : synthesis-to-sequencing cost ratio
        beta     : strand-length scaling factor (must be > 1)
        c_min    : lower bound of search grid
        c_max    : upper bound of search grid
        n_points : number of grid points

    Returns:
        c* as a float, or math.nan if beta <= 1 (no positive rate achievable).
    """
    if beta <= 1:
        return math.nan

    c_grid = np.linspace(c_min, c_max, n_points)
    costs  = (q + c_grid) / (1 - np.exp(-c_grid))   # (1-1/beta) cancels in argmin
    return float(c_grid[np.argmin(costs)])
