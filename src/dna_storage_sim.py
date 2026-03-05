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
import csv
import math
import os
import numpy as np
import matplotlib.pyplot as plt

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

    Q = seen.sum(axis=1)              # unique strands recovered per trial
    Q_over_M    = Q / M               # fraction recovered, one value per trial
    Rs          = Q_over_M * (1 - 1 / beta)   # storage rate per unique strand
    Rs_physical = Rs / DUP            # storage rate per physical strand

    n = trials
    mean_Q  = float(Q_over_M.mean())
    se_Q    = float(Q_over_M.std(ddof=1) / math.sqrt(n)) if n > 1 else 0.0
    mean_Rs = float(Rs.mean())
    se_Rs   = float(Rs.std(ddof=1)   / math.sqrt(n)) if n > 1 else 0.0
    mean_Rp = float(Rs_physical.mean())

    return {
        "M": M,
        "L": L,
        "N": N,
        "total_physical": total_physical,
        "mean_Q_over_M": mean_Q,
        "se_Q_over_M": se_Q,
        # kept for backwards-compat with print code
        "avg_fraction": mean_Q,
        "std_fraction": float(Q_over_M.std(ddof=1)) if n > 1 else 0.0,
        "theoretical": theoretical,
        "mean_Rs": mean_Rs,
        "se_Rs": se_Rs,
        "mean_Rs_physical": mean_Rp,
        "mode": mode,
    }


def plot_results(csv_rows: list, limit: float, beta: float, mode: str, c: float,
                 plots_dir: str = "../plots"):
    """
    Generate and save two plots from sweep results.

    Plot 1: Fraction recovered vs M (log x-axis)
    Plot 2: Storage rates (Rs, Rs_physical) vs M
    """
    os.makedirs(plots_dir, exist_ok=True)

    Ms            = [r["M"]                for r in csv_rows]
    mean_Q        = [r["mean_Q_over_M"]    for r in csv_rows]
    se_Q          = [r["se_Q_over_M"]      for r in csv_rows]
    mean_Rs       = [r["mean_Rs"]          for r in csv_rows]
    se_Rs         = [r["se_Rs"]            for r in csv_rows]
    mean_Rp       = [r["mean_Rs_physical"] for r in csv_rows]
    theory_Rs     = [limit * (1 - 1 / beta)] * len(Ms)

    # --- Plot 1: Fraction recovered vs M ---
    fig, ax = plt.subplots()
    ax.errorbar(Ms, mean_Q, yerr=se_Q, fmt="o-", capsize=4, label="Simulated")
    ax.axhline(limit, color="red", linestyle="--", label=f"Theory = {limit:.4f}")
    ax.set_xscale("log")
    ax.set_xlabel("M (number of unique strands)")
    ax.set_ylabel("Fraction of unique strands recovered")
    ax.set_title(f"Fraction Recovered vs M  [mode={mode}, c={c}, DUP={DUP}]")
    ax.legend()
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "fraction_recovered.png"), dpi=150)
    plt.close(fig)
    print(f"Saved {plots_dir}/fraction_recovered.png")

    # --- Plot 2: Storage rates vs M ---
    fig, ax = plt.subplots()
    ax.errorbar(Ms, mean_Rs, yerr=se_Rs, fmt="s-", capsize=4, label="Rs (per unique)")
    ax.plot(Ms, mean_Rp, "^--", label="Rs_physical (per physical)")
    ax.plot(Ms, theory_Rs, "r:", label=f"Theory Rs = {theory_Rs[0]:.4f}")
    ax.set_xscale("log")
    ax.set_xlabel("M (number of unique strands)")
    ax.set_ylabel("Storage rate")
    ax.set_title(f"Storage Rate vs M  [mode={mode}, c={c}, beta={beta}, DUP={DUP}]")
    ax.legend()
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "storage_rate.png"), dpi=150)
    plt.close(fig)
    print(f"Saved {plots_dir}/storage_rate.png")


def main():
    parser = argparse.ArgumentParser(
        description="DNA Storage Recovery Simulator",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--M_values", type=int, nargs="+",
        default=[
            10, 13, 18, 25, 35, 48, 67, 92, 126, 174,
            239, 329, 452, 621, 853, 1172, 1610, 2212, 3039, 4175,
            5736, 7880, 10826, 14873, 20433, 28072, 38566, 52983, 72789, 100000,
        ],
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
    parser.add_argument(
        "--output", type=str, default=None,
        help="Path to CSV file for saving sweep results (optional)",
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

    col = "{:>8}  {:>6}  {:>10}  {:>10}  {:>11}  {:>8}  {:>8}  {:>10}  {:>11}"
    print(col.format(
        "M", "L", "N reads",
        "Simulated", "Theoretical", "Error", "Std",
        "Rs", "Rs_physical",
    ))
    print(col.format(
        "-"*8, "-"*6, "-"*10,
        "-"*10, "-"*11, "-"*8, "-"*8,
        "-"*10, "-"*11,
    ))

    CSV_FIELDS = [
        "M", "c", "beta", "DUP", "mode",
        "mean_Q_over_M", "se_Q_over_M",
        "mean_Rs", "se_Rs",
        "mean_Rs_physical",
        "theory_fraction",
    ]

    csv_rows = []

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
            f"{r['mean_Rs']:.6f}",
            f"{r['mean_Rs_physical']:.6f}",
        ))

        csv_rows.append({
            "M": M,
            "c": args.c,
            "beta": args.beta,
            "DUP": DUP,
            "mode": args.mode,
            "mean_Q_over_M": r["mean_Q_over_M"],
            "se_Q_over_M": r["se_Q_over_M"],
            "mean_Rs": r["mean_Rs"],
            "se_Rs": r["se_Rs"],
            "mean_Rs_physical": r["mean_Rs_physical"],
            "theory_fraction": r["theoretical"],
        })

    print()
    print(f"As M grows, simulated fraction converges to {limit:.6f}")

    if args.output:
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CSV_FIELDS)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"\nResults saved to {args.output}")

    plot_results(csv_rows, limit, args.beta, args.mode, args.c,
                 plots_dir=os.path.join(os.path.dirname(__file__), "../plots"))


if __name__ == "__main__":
    main()
