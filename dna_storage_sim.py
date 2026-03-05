"""
DNA Storage Recovery Simulator

Simulates recovery of unique DNA strands during sequencing and verifies
that the fraction recovered approaches a theoretical limit as M grows.

Model:
  - M unique DNA strands, each duplicated DUP times before sequencing
  - Reads sample strands uniformly at random from the physical pool
  - We measure how many unique strands are recovered

Modes:
  per_unique  : N = c * M reads (coverage relative to unique strands)
  per_physical: N = c * M * DUP reads (coverage relative to physical pool)

Theoretical expected fraction recovered:
  P(a given unique strand is missed) = (1 - 1/M)^N  ~  e^(-N/M)
  E[fraction recovered]              = 1 - e^(-N/M)

  per_unique   -> 1 - e^(-c)
  per_physical -> 1 - e^(-c * DUP)
"""

import argparse
import math
import numpy as np

DUP = 10


def simulate_recovery(M: int, c: float, beta: float, trials: int, mode: str) -> dict:
    """
    Run `trials` independent recovery simulations for a given M.

    Args:
        M      : Number of unique strands
        c      : Coverage constant
        beta   : Strand-length scaling factor (L = beta * log2(M) nucleotides)
        trials : Number of independent trials
        mode   : "per_unique" or "per_physical"

    Returns:
        dict with keys: M, L, N, total_physical, avg_fraction, std_fraction,
                        theoretical, mode
    """
    L = max(1, int(beta * math.log2(M)))   # strand length (nucleotides)
    total_physical = M * DUP

    if mode == "per_unique":
        N = round(c * M)
        theoretical = 1 - math.exp(-c)
    elif mode == "per_physical":
        N = round(c * total_physical)
        theoretical = 1 - math.exp(-c * DUP)
    else:
        raise ValueError(f"Unknown mode '{mode}'. Choose 'per_unique' or 'per_physical'.")

    rng = np.random.default_rng()

    # Sample all trials at once: shape (trials, N), values in [0, M)
    # No list of size DUP*M is ever created; equal multiplicity means
    # sampling a physical strand is identical to sampling a unique strand ID.
    reads = rng.integers(0, M, size=(trials, N))

    # Boolean presence array: True where strand ID appears in each trial
    seen = np.zeros((trials, M), dtype=bool)
    seen[np.arange(trials)[:, None], reads] = True

    Q = seen.sum(axis=1)          # unique strands recovered per trial
    Q_over_M = Q / M              # fraction recovered per trial

    return {
        "M": M,
        "L": L,
        "N": N,
        "total_physical": total_physical,
        "avg_fraction": float(Q_over_M.mean()),
        "std_fraction": float(Q_over_M.std(ddof=1)) if trials > 1 else 0.0,
        "theoretical": theoretical,
        "mode": mode,
    }


def main():
    parser = argparse.ArgumentParser(
        description="DNA Storage Recovery Simulator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--M_values", type=int, nargs="+", default=[10, 50, 100, 500, 1000, 5000],
        help="List of M values (number of unique strands)",
    )
    parser.add_argument(
        "--c", type=float, default=1.0,
        help="Coverage constant",
    )
    parser.add_argument(
        "--beta", type=float, default=1.0,
        help="Strand-length scaling factor: L = beta * log2(M) nucleotides",
    )
    parser.add_argument(
        "--trials", type=int, default=50,
        help="Number of independent simulation trials per M value",
    )
    parser.add_argument(
        "--mode", type=str, choices=["per_unique", "per_physical"],
        default="per_unique",
        help=(
            "Sampling mode: "
            "per_unique => N = c*M reads; "
            "per_physical => N = c*M*DUP reads"
        ),
    )
    args = parser.parse_args()

    # Theoretical limit as M -> inf
    if args.mode == "per_unique":
        limit = 1 - math.exp(-args.c)
        limit_expr = f"1 - e^(-c) = 1 - e^(-{args.c})"
    else:
        limit = 1 - math.exp(-args.c * DUP)
        limit_expr = f"1 - e^(-c*DUP) = 1 - e^(-{args.c}*{DUP})"

    # Header
    print()
    print("DNA Storage Recovery Simulator")
    print("=" * 70)
    print(f"  DUP    = {DUP}  (copies per unique strand)")
    print(f"  c      = {args.c}")
    print(f"  beta   = {args.beta}  (L = beta * log2(M) nucleotides)")
    print(f"  mode   = {args.mode}")
    print(f"  trials = {args.trials}")
    print(f"  Theoretical limit (M->inf): {limit_expr} = {limit:.6f}")
    print("=" * 70)

    col = "{:>8}  {:>6}  {:>10}  {:>12}  {:>12}  {:>8}  {:>8}"
    print(col.format("M", "L", "N reads", "Simulated", "Theoretical", "Error", "Std"))
    print(col.format("-" * 8, "-" * 6, "-" * 10, "-" * 12, "-" * 12, "-" * 8, "-" * 8))

    for M in args.M_values:
        r = simulate_recovery(M, args.c, args.beta, args.trials, args.mode)
        error = abs(r["avg_fraction"] - r["theoretical"])
        print(col.format(
            f"{M:,}",
            r["L"],
            f"{r['N']:,}",
            f"{r['avg_fraction']:.6f}",
            f"{r['theoretical']:.6f}",
            f"{error:.6f}",
            f"{r['std_fraction']:.6f}",
        ))

    print()
    print(f"As M grows, simulated fraction converges to {limit:.6f}")


if __name__ == "__main__":
    main()
