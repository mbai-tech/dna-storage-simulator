import argparse
import csv
import math
import os
import numpy as np
import matplotlib.pyplot as plt

from cost_model import Rs_theory, Rr_theory, cost_theory

DUP = 10


def simulate_recovery(M: int, c: float, beta: float, trials: int, mode: str,
                      target_recovery: float = None) -> dict:
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

    # --- Runtime and memory ---
    # Time per trial: O(N) where N = c*M (per_unique) or c*M*DUP (per_physical).
    #   - Drawing N integers costs O(N).
    #   - Marking the seen array via scatter-index is also O(N).
    #   - There is no loop over molecules, no list of size DUP*M, and no
    #     Python-level iteration over reads — all trials are vectorised into
    #     one numpy call of shape (trials, N).
    # Runtime therefore scales linearly with the number of reads N, and
    # since N = c*M the dominant cost is simply proportional to M.
    # This matches the design goal: runtime depends mostly on M and grows
    # linearly with coverage c.
    #
    # Memory per trial: O(M) for the boolean "seen" array of length M.
    #   - The reads array is (trials, N) = O(trials * c * M), but is freed
    #     immediately after the scatter step.
    #   - The seen array is (trials, M) booleans = trials * M bytes.
    reads = rng.integers(0, M, size=(trials, N))

    # Boolean presence array: True where strand ID appears in each trial.
    # Shape (trials, M) — one bool per unique strand per trial.
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

    # Reliability: lower confidence bound must meet target_recovery
    reliable = None
    if target_recovery is not None:
        reliable = (mean_Q - 2 * se_Q) >= target_recovery

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
        "reliable": reliable,
        "mode": mode,
    }


def find_c_star(q: float, M: int, beta: float, target_recovery: float,
                trials: int, mode: str = "per_unique") -> dict:
    """
    Find the coverage c* that minimises cost(c, q, beta) subject to the
    reliability constraint  mean_Q - 2*se_Q >= target_recovery.

    Search strategy:
      Phase 1 — coarse grid  c in [0.2, 20.0] step 0.1
      Phase 2 — fine grid    c in [c_coarse - 0.1, c_coarse + 0.1] step 0.01

    Args:
        q               : synthesis-to-sequencing cost ratio
        M               : number of unique strands
        beta            : strand-length scaling factor
        target_recovery : minimum acceptable mean_Q - 2*se_Q
        trials          : simulation trials per candidate c
        mode            : "per_unique" or "per_physical"

    Returns:
        dict with keys: c_star, best_cost, mean_Q_over_M, se_Q_over_M.
        c_star is None if no reliable c was found on the grid.
    """
    def _evaluate(c_vals):
        """Return (c_star, best_cost, best_r) over an iterable of c values."""
        best_c    = None
        best_cost = math.inf
        best_r    = None
        for c in c_vals:
            r = simulate_recovery(M, float(c), beta, trials, mode,
                                  target_recovery=target_recovery)
            if r["reliable"]:
                cst = cost_theory(float(c), q, beta)
                if cst < best_cost:
                    best_cost = cst
                    best_c    = float(c)
                    best_r    = r
        return best_c, best_cost, best_r

    # Phase 1: coarse grid
    coarse_grid = np.arange(0.2, 20.1, 0.1)
    best_c, best_cost, best_r = _evaluate(coarse_grid)

    # Phase 2: refinement around coarse best
    if best_c is not None:
        fine_grid = np.arange(max(0.1, best_c - 0.1),
                              min(20.0, best_c + 0.1) + 0.001,
                              0.01)
        refined_c, refined_cost, refined_r = _evaluate(fine_grid)
        if refined_c is not None and refined_cost < best_cost:
            best_c, best_cost, best_r = refined_c, refined_cost, refined_r

    if best_r is None:
        return {"c_star": None, "best_cost": math.inf,
                "mean_Q_over_M": None, "se_Q_over_M": None}

    return {
        "c_star":        best_c,
        "best_cost":     best_cost,
        "mean_Q_over_M": best_r["mean_Q_over_M"],
        "se_Q_over_M":   best_r["se_Q_over_M"],
    }


def plot_results(csv_rows: list, limit: float, beta: float, mode: str, c: float,
                 plots_dir: str = "../plots", target_recovery: float = None):
    """
    Generate and save two plots from sweep results.

    Plot 1: Fraction recovered vs M (log x-axis)
    Plot 2: Storage rates (Rs, Rs_physical) vs M
    """
    os.makedirs(plots_dir, exist_ok=True)

    Ms            = [r["M"]                for r in csv_rows]
    mean_Q        = [r["mean_Q_over_M"]    for r in csv_rows]
    se_Q          = [r["se_Q_over_M"]      for r in csv_rows]
    mean_Rp       = [r["mean_Rs_physical"] for r in csv_rows]

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

    # --- Plot 3: Storage rate (per physical strand) vs M ---
    fig, ax = plt.subplots()
    ax.plot(Ms, mean_Rp, "o-", color="tab:blue", label="Rs_physical (per physical)")
    theory_Rp = [limit * (1 - 1 / beta) / DUP] * len(Ms)
    ax.plot(Ms, theory_Rp, "--", color="red", label=f"Theory Rs_physical = {theory_Rp[0]:.4f}")
    ax.set_xscale("log")
    ax.set_xlabel("M (number of unique strands)")
    ax.set_ylabel("Storage rate (bits per nucleotide synthesized)")
    ax.set_title(f"Storage Rate (Physical) vs M  [mode={mode}, c={c}, beta={beta}, DUP={DUP}]")
    ax.legend()
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "storage_rate_physical.png"), dpi=150)
    plt.close(fig)
    print(f"Saved {plots_dir}/storage_rate_physical.png")


def plot_design_recommendations(q_rows: list, target_recovery: float,
                                plots_dir: str = "../plots"):
    """
    Generate and save two design-recommendation plots from a q-sweep.

    Plot 1: c* vs q  (one line per unique M value, log x-axis)
    Plot 2: best_cost vs q  (one line per unique M value, log x-axis)

    Args:
        q_rows         : list of row dicts produced by the q-sweep loop
        target_recovery: reliability target used in the sweep
        plots_dir      : directory to save PNG files
    """
    os.makedirs(plots_dir, exist_ok=True)

    # Group rows by M
    M_vals = sorted(set(r["M"] for r in q_rows))
    q_vals = sorted(set(r["q"] for r in q_rows))

    def _series(M, key):
        subset = [r for r in q_rows if r["M"] == M]
        subset.sort(key=lambda r: r["q"])
        qs  = [r["q"]  for r in subset if r[key] is not None]
        ys  = [r[key]  for r in subset if r[key] is not None]
        return qs, ys

    # --- Plot 1: c* vs q ---
    fig, ax = plt.subplots()
    for M in M_vals:
        qs, cs = _series(M, "c_star")
        ax.plot(qs, cs, "o-", label=f"M={M:,}")
    ax.set_xscale("log")
    ax.set_xlabel("q  (synthesis / sequencing cost ratio)")
    ax.set_ylabel("Optimal coverage c*")
    ax.set_title(
        f"Design Recommendation: c* vs q\n"
        f"[target_recovery={target_recovery}, DUP={DUP}]"
    )
    ax.legend(fontsize=8)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    fig.tight_layout()
    path1 = os.path.join(plots_dir, "design_c_star_vs_q.png")
    fig.savefig(path1, dpi=150)
    plt.close(fig)
    print(f"Saved {path1}")

    # --- Plot 2: best_cost vs q ---
    fig, ax = plt.subplots()
    for M in M_vals:
        qs, costs = _series(M, "best_cost")
        ax.plot(qs, costs, "s-", label=f"M={M:,}")
    ax.set_xscale("log")
    ax.set_xlabel("q  (synthesis / sequencing cost ratio)")
    ax.set_ylabel("Minimum cost per stored bit  (q/Rs + 1/Rr)")
    ax.set_title(
        f"Design Recommendation: Cost* vs q\n"
        f"[target_recovery={target_recovery}, DUP={DUP}]"
    )
    ax.legend(fontsize=8)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    fig.tight_layout()
    path2 = os.path.join(plots_dir, "design_cost_vs_q.png")
    fig.savefig(path2, dpi=150)
    plt.close(fig)
    print(f"Saved {path2}")


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
    parser.add_argument(
        "--q", type=float, default=10000.0,
        help="Synthesis-to-sequencing cost ratio (paper uses q ~ 10,000-100,000)",
    )
    parser.add_argument(
        "--target_recovery", type=float, default=None,
        help=(
            "Reliability target: (c, M) is reliable if "
            "mean_Q_over_M - 2*se_Q_over_M >= target_recovery. "
            "Omit to skip reliability check."
        ),
    )
    parser.add_argument(
        "--find_c_star", action="store_true",
        help=(
            "Search over c in [0.2, 20] to find the coverage that minimises "
            "cost(c, q, beta) subject to the reliability constraint. "
            "Requires --target_recovery."
        ),
    )
    parser.add_argument(
        "--q_values", type=float, nargs="+", default=None,
        help=(
            "List of synthesis-to-sequencing cost ratios to sweep. "
            "For each q and each M, runs find_c_star and saves one CSV row. "
            "Paper notes typical q ~ 10,000-100,000. Requires --target_recovery."
        ),
    )
    parser.add_argument(
        "--q_output", type=str, default="q_sweep.csv",
        help="CSV file path for --q_values sweep results (default: q_sweep.csv).",
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
    print()
    print("  Cost model (Heckel et al. 2017, arXiv:1705.04732):")
    print(f"    q (synthesis/sequencing cost ratio) = {args.q}")
    print(f"    Rs_theory  = (1 - e^-c)(1 - 1/beta) = {Rs_theory(args.c, args.beta):.6f}")
    print(f"    Rr_theory  = Rs_theory / c           = {Rr_theory(args.c, args.beta):.6f}")
    print(f"    cost/bit   = q/Rs + 1/Rr             = {cost_theory(args.c, args.q, args.beta):.4f}")
    if args.target_recovery is not None:
        print(f"  Reliability target: mean - 2*SE >= {args.target_recovery}")
    print("=" * 70)

    if args.find_c_star:
        if args.target_recovery is None:
            print("WARNING: --find_c_star requires --target_recovery; skipping.")
        else:
            print()
            print("Optimal coverage search (find_c_star)")
            print("-" * 70)
            cstar_col = "{:>8}  {:>8}  {:>12}  {:>12}  {:>10}"
            print(cstar_col.format("M", "c*", "mean_Q", "se_Q", "cost*"))
            print(cstar_col.format("-"*8, "-"*8, "-"*12, "-"*12, "-"*10))
            for M in args.M_values:
                res = find_c_star(args.q, M, args.beta, args.target_recovery,
                                  args.trials, args.mode)
                if res["c_star"] is None:
                    print(cstar_col.format(f"{M:,}", "none", "—", "—", "inf"))
                else:
                    print(cstar_col.format(
                        f"{M:,}",
                        f"{res['c_star']:.2f}",
                        f"{res['mean_Q_over_M']:.6f}",
                        f"{res['se_Q_over_M']:.6f}",
                        f"{res['best_cost']:.2f}",
                    ))
            print()

    if args.q_values:
        if args.target_recovery is None:
            print("WARNING: --q_values requires --target_recovery; skipping q sweep.")
        else:
            Q_SWEEP_FIELDS = [
                "q", "M", "beta", "DUP", "target_recovery",
                "c_star", "best_cost",
                "mean_Q_over_M_at_c_star", "se_Q_over_M_at_c_star",
                "Rs_theory_at_c_star", "Rr_theory_at_c_star",
            ]
            q_rows = []

            print()
            print("q-sweep (find_c_star for each q × M)")
            print("-" * 80)
            qcol = "{:>10}  {:>8}  {:>8}  {:>12}  {:>12}  {:>14}  {:>14}"
            print(qcol.format("q", "M", "c*", "mean_Q", "se_Q", "Rs_theory", "cost*"))
            print(qcol.format("-"*10, "-"*8, "-"*8, "-"*12, "-"*12, "-"*14, "-"*14))

            for q in args.q_values:
                for M in args.M_values:
                    res = find_c_star(q, M, args.beta, args.target_recovery,
                                      args.trials, args.mode)
                    c_s = res["c_star"]
                    rs  = Rs_theory(c_s, args.beta) if c_s is not None else None
                    rr  = Rr_theory(c_s, args.beta) if c_s is not None else None

                    print(qcol.format(
                        f"{q:,.0f}",
                        f"{M:,}",
                        f"{c_s:.2f}"   if c_s is not None else "none",
                        f"{res['mean_Q_over_M']:.6f}" if c_s is not None else "—",
                        f"{res['se_Q_over_M']:.6f}"   if c_s is not None else "—",
                        f"{rs:.6f}"    if rs  is not None else "—",
                        f"{res['best_cost']:.2f}"     if c_s is not None else "inf",
                    ))

                    q_rows.append({
                        "q":                       q,
                        "M":                       M,
                        "beta":                    args.beta,
                        "DUP":                     DUP,
                        "target_recovery":         args.target_recovery,
                        "c_star":                  c_s,
                        "best_cost":               res["best_cost"],
                        "mean_Q_over_M_at_c_star": res["mean_Q_over_M"],
                        "se_Q_over_M_at_c_star":   res["se_Q_over_M"],
                        "Rs_theory_at_c_star":     rs,
                        "Rr_theory_at_c_star":     rr,
                    })

            with open(args.q_output, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=Q_SWEEP_FIELDS)
                writer.writeheader()
                writer.writerows(q_rows)
            print(f"\nq-sweep results saved to {args.q_output}")

            # Human-readable design summary
            print()
            print("Design Recommendations")
            print("-" * 70)
            for r in q_rows:
                c_s = r["c_star"]
                if c_s is None:
                    print(
                        f"  q={r['q']:>10,.0f}, M={r['M']:>8,}: "
                        f"no reliable c found (target {args.target_recovery} unreachable)"
                    )
                else:
                    met = r["mean_Q_over_M_at_c_star"] is not None
                    print(
                        f"  q={r['q']:>10,.0f}, M={r['M']:>8,}: "
                        f"optimal c* = {c_s:.2f}, "
                        f"recovery = {r['mean_Q_over_M_at_c_star']:.4f} "
                        f"(target {args.target_recovery}), "
                        f"cost = {r['best_cost']:.1f}"
                    )
            print()

            plots_dir = os.path.join(os.path.dirname(__file__), "../plots")
            plot_design_recommendations(q_rows, args.target_recovery,
                                        plots_dir=plots_dir)
            print()

    col = "{:>8}  {:>6}  {:>10}  {:>10}  {:>11}  {:>8}  {:>8}  {:>10}  {:>11}  {:>10}"
    print(col.format(
        "M", "L", "N reads",
        "Simulated", "Theoretical", "Error", "Std",
        "Rs", "Rs_physical", "Reliable",
    ))
    print(col.format(
        "-"*8, "-"*6, "-"*10,
        "-"*10, "-"*11, "-"*8, "-"*8,
        "-"*10, "-"*11, "-"*10,
    ))

    CSV_FIELDS = [
        "M", "c", "beta", "DUP", "mode",
        "mean_Q_over_M", "se_Q_over_M",
        "mean_Rs", "se_Rs",
        "mean_Rs_physical",
        "theory_fraction",
        "Rs_theory", "Rr_theory", "cost_theory",
        "reliable",
    ]

    csv_rows = []

    for M in args.M_values:
        r = simulate_recovery(M, args.c, args.beta, args.trials, args.mode,
                              target_recovery=args.target_recovery)
        error = abs(r["avg_fraction"] - r["theoretical"])

        if r["reliable"] is None:
            rel_str = "N/A"
        else:
            rel_str = "YES" if r["reliable"] else "NO"

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
            rel_str,
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
            "Rs_theory": Rs_theory(args.c, args.beta),
            "Rr_theory": Rr_theory(args.c, args.beta),
            "cost_theory": cost_theory(args.c, args.q, args.beta),
            "reliable": r["reliable"],
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
                 plots_dir=os.path.join(os.path.dirname(__file__), "../plots"),
                 target_recovery=args.target_recovery)


if __name__ == "__main__":
    main()
