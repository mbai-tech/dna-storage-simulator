# DNA Storage Simulator

A Python simulation of DNA-based data storage, modeling sequencing recovery rates and verifying them against theoretical limits from [Heckel et al. 2017 (arXiv:1705.04732)](https://arxiv.org/abs/1705.04732).

---

## Setup

```bash
python3 -m venv venv
source venv/bin/activate
pip install numpy matplotlib
```

## Usage

```bash
source venv/bin/activate

# Basic run
python src/dna_storage_sim.py --beta 2.0

# With reliability check
python src/dna_storage_sim.py --beta 2.0 --target_recovery 0.60

# Find optimal c* for a given cost ratio
python src/dna_storage_sim.py --beta 2.0 --target_recovery 0.60 --find_c_star --q 10000

# q-sweep with design recommendation plots
python src/dna_storage_sim.py --beta 2.0 --target_recovery 0.60 \
  --q_values 1 100 1000 10000 100000 --q_output q_sweep.csv
```

## CLI Options

| Flag | Default | Description |
|---|---|---|
| `--M_values` | 30 log-spaced points (10–100,000) | Number of unique strands |
| `--c` | 1.0 | Coverage depth (reads per unique strand) |
| `--beta` | 1.0 | Strand-length factor: L = beta × log₂(M); use > 1 |
| `--trials` | 50 | Simulation trials per M value |
| `--mode` | `per_unique` | `per_unique` (N=c×M) or `per_physical` (N=c×M×DUP) |
| `--q` | 10000 | Synthesis-to-sequencing cost ratio |
| `--target_recovery` | None | Reliability threshold: mean − 2×SE ≥ target |
| `--find_c_star` | off | Find cost-minimising coverage c* |
| `--q_values` | None | Sweep multiple q values; generates design plots |
| `--q_output` | `q_sweep.csv` | Output CSV for q-sweep |
| `--output` | None | Output CSV for main sweep |

## Model

M unique strands, each duplicated DUP=10 times. N reads sample the pool uniformly at random. Recovery fraction per trial = (unique strands seen) / M.

**Theoretical expected fraction recovered:**

| Mode | N | Theory |
|---|---|---|
| `per_unique` | c × M | 1 − e^{−c} |
| `per_physical` | c × M × DUP | 1 − e^{−c × DUP} |

## Cost Model

Based on Heckel et al. 2017, Section 4:

| Symbol | Formula | Meaning |
|---|---|---|
| `Rs` | `(1 − e^{−c})(1 − 1/β)` | Bits stored per nucleotide synthesized |
| `Rr` | `Rs / c` | Bits recovered per nucleotide sequenced |
| `cost` | `q/Rs + 1/Rr` | Cost per stored bit |

For q = 10,000, the optimal coverage minimizing cost is **c* ≈ 9.2**.

