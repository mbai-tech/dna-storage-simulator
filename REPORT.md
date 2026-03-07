# DNA-Based Data Storage: Simulation and Cost Analysis

**Bi 001C Final Project**

---

## 1. Introduction

DNA is an increasingly attractive medium for long-term data storage. It offers extraordinary information density (theoretically up to ~2 bits per nucleotide), physical stability over thousands of years, and an exponentially declining synthesis cost. However, retrieving data from a DNA pool requires sequencing, which introduces a fundamental probabilistic problem: not every molecule in the pool will be read. Some unique sequences may be missed entirely, leading to data loss.

This project implements a simulation of DNA-based data storage with two goals:

1. **Verify** that the fraction of unique strands recovered during sequencing matches the theoretical prediction from the coupon-collector model.
2. **Evaluate** the cost tradeoffs between synthesis and sequencing using the framework from Heckel et al. (2017), and find the coverage depth that minimizes total cost subject to a reliability constraint.

The theoretical framework is drawn from:

> R. Heckel, G. Mikutis, and R. N. Grass, "A Characterization of the DNA Data Storage Channel," *arXiv:1705.04732*, 2017.

---

## 2. Model

### 2.1 Physical Setup

The storage system is modeled as follows:

- **M** unique DNA strands encode the data, each L nucleotides long, where L = β · log₂(M).
- Each unique strand is physically duplicated **DUP = 10** times before sequencing, yielding a pool of M × DUP total molecules.
- During sequencing, **N** reads are drawn uniformly at random from the physical pool. Because all DUP copies of a strand are identical, sampling any copy is equivalent to sampling the unique strand ID. Under the `per_unique` mode, the duplication cancels out and the effective pool is just the M unique IDs.

### 2.2 Coverage Modes

| Mode | N (number of reads) | Effective pool |
|---|---|---|
| `per_unique` | c × M | M unique strands |
| `per_physical` | c × M × DUP | M × DUP physical strands |

The parameter **c** is the coverage depth: the average number of times each unique strand is expected to be read.

### 2.3 Theoretical Recovery Fraction

This is a coupon-collector problem. Each read samples a unique strand ID uniformly in [0, M). The probability that a given strand is missed in N reads is:

```
P(missed) = (1 − 1/M)^N  ≈  e^{−N/M}
```

So the expected fraction of strands recovered is:

```
E[fraction recovered] = 1 − e^{−N/M}
```

Substituting the two modes:

| Mode | Theory |
|---|---|
| `per_unique` | 1 − e^{−c} |
| `per_physical` | 1 − e^{−c · DUP} |

As M grows large, the simulated recovery fraction converges to this theoretical limit.

---

## 3. Simulation

### 3.1 Implementation

The simulator (`src/dna_storage_sim.py`) runs `trials` independent recovery experiments for each value of M. Each trial:

1. Draws N strand IDs uniformly at random using `numpy.random.integers`.
2. Records which unique strands were seen using a boolean array of length M.
3. Computes Q/M, the fraction of unique strands recovered.

All trials are vectorized into a single NumPy operation of shape `(trials, N)`, avoiding Python-level loops. Runtime scales as O(N) = O(c · M) per trial; memory is O(trials · M).

### 3.2 Verified Convergence

Running the simulator with β = 2.0, c = 1.0, and 50 trials across 30 log-spaced values of M from 10 to 100,000, the simulated recovery fraction converges to the theoretical limit of 1 − e⁻¹ ≈ 0.6321 as M grows. At small M, the coupon-collector approximation is less accurate due to finite-size effects, but by M ~ 10,000 the simulation matches theory to within ±0.001.

### 3.3 Reliability Rule

For a given (c, M) configuration to be considered **reliable**, its simulated recovery must satisfy a lower confidence bound:

```
mean(Q/M) − 2 · SE(Q/M)  ≥  target_recovery
```

where SE is the standard error across trials. This is a conservative 2-sigma rule: if the lower bound of the confidence interval exceeds the target, we are confident the system reliably meets it.

---

## 4. Cost Model

### 4.1 Rates and Cost (Heckel et al. 2017)

Heckel et al. define two efficiency measures and a cost function (Section 4, Theorem 1, Corollary 1):

| Symbol | Formula | Meaning |
|---|---|---|
| Rs | (1 − e^{−c})(1 − 1/β) | Bits stored per nucleotide **synthesized** |
| Rr | Rs / c | Bits recovered per nucleotide **sequenced** |
| cost | q/Rs + 1/Rr | Total cost per stored bit |

The factor (1 − 1/β) accounts for the overhead of using strand length L = β · log₂(M) to index M strands — only (1 − 1/β) of each strand's bits carry payload data.

The synthesis-to-sequencing cost ratio **q** reflects the fact that, in practice, DNA synthesis is far more expensive than sequencing. Current estimates place q in the range of 10,000–100,000.

### 4.2 Optimal Coverage c*

Differentiating cost with respect to c yields the optimality condition:

```
e^c = q + c + 1
```

This has no closed-form solution but a unique numerical solution for any q > 0. The optimal c* depends only on q, not on β or M. For q = 10,000, this gives **c* ≈ 9.2**, matching the paper's result in Section 4.

The simulator finds c* via a two-phase grid search:
- Phase 1: coarse grid over c ∈ [0.2, 20.0] at step 0.1
- Phase 2: fine grid over ±0.1 around the coarse optimum at step 0.01

The search is constrained to c values where the reliability rule is satisfied.

---

## 5. Results

### 5.1 Recovery vs. M

The fraction-recovered plot (β = 2.0, c = 2.0, 50 trials) shows the simulated recovery converging to the theoretical limit as M grows. At M = 10, stochastic variance is high. By M = 10,000, the simulated mean is within 0.1% of theory. This confirms that the large-M approximation underlying the Heckel et al. model is valid.

### 5.2 Optimal Coverage vs. Cost Ratio

The q-sweep was run with β = 2.0, target recovery 0.60, across q ∈ {1, 100, 500, 1000, 2500, 5000, 10000, 25000, 50000, 75000, 100000} and M ∈ {1000, 10000}. Selected results:

| q | c* | Rs at c* | Cost per bit |
|---|---|---|---|
| 1 | 1.15 | 0.342 | 6.29 |
| 100 | 4.66 | 0.495 | 211 |
| 1,000 | 6.92 | 0.500 | 2,016 |
| 10,000 | 9.21 | 0.500 | 20,020 |
| 100,000 | 11.51 | 0.500 | 200,025 |

Several observations:

- **c* increases with q.** When synthesis is expensive relative to sequencing, it is worth sequencing more (higher c) to maximize the return on synthesized strands. The relationship is approximately logarithmic.
- **Rs saturates near 0.5 for β = 2.0.** The theoretical maximum is (1 − 1/β) = 0.5, achieved as c → ∞. By c* ≈ 9, recovery is > 99.9% and Rs is within 0.1% of its ceiling.
- **Cost scales roughly linearly with q.** Since Rs saturates, the dominant term in cost = q/Rs + 1/Rr grows linearly with q at large q.
- **Results are consistent across M = 1,000 and M = 10,000**, confirming that c* depends on q but not on M, as predicted by the theory.

### 5.3 Recovery at c*

At q = 10,000 and M = 10,000, the simulator finds c* = 9.21 with a simulated recovery of 99.99% (mean Q/M = 0.99988, SE = 0.000019), far exceeding the 60% reliability target. This illustrates a key point: the optimal c is not the minimum c needed for reliability — it is the c that minimizes total cost given the constraint. Because Rs saturates, the marginal cost of extra coverage grows faster than its benefit, making very high coverage wasteful.

---

## 6. Discussion

### 6.1 The Coverage Tradeoff

The cost model reveals a fundamental tradeoff: low coverage wastes synthesized molecules (many unique strands go unread), while high coverage wastes sequencing reads (diminishing returns as nearly all strands have been seen). The optimal c* balances these two costs. At q = 10,000, c* ≈ 9.2 represents approximately 9 reads per unique strand — notably higher than the minimum coverage needed to achieve reliable recovery.

### 6.2 The Role of β

The strand-length factor β controls information density. With β = 1, each strand's entire length is used for addressing — no payload capacity, and Rs = 0. As β increases, more of each strand carries data, and Rs approaches 1 − e^{-c}. In practice β must be tuned to the desired library size M: β = 2.0 means strands are twice as long as the minimum needed to address M unique sequences.

### 6.3 Limitations

- **Uniform sampling.** The model assumes reads are drawn uniformly from the pool. In practice, PCR amplification and sequencing introduce non-uniformity and coverage bias.
- **No error model.** Sequencing errors (substitutions, insertions, deletions) are not simulated. Heckel et al. treat these as a separate channel impairment.
- **Fixed DUP.** The number of physical copies per strand is fixed at 10. In practice this is controlled by PCR cycle count and affects the actual pool distribution.
- **Large-M approximation.** The theoretical formulas are asymptotic in M. The simulation shows they hold well by M ~ 10,000, but may not apply to very small libraries.

---

## 7. Conclusion

This project demonstrates, through direct Monte Carlo simulation, that the fraction of unique DNA strands recovered during sequencing follows the coupon-collector formula 1 − e^{−c}, converging to the theoretical limit as M grows. Using the cost model of Heckel et al. (2017), we find that the optimal sequencing coverage c* increases logarithmically with the synthesis-to-sequencing cost ratio q, reaching c* ≈ 9.2 at the practically relevant value of q = 10,000. These results confirm the paper's key theoretical predictions and provide a practical framework for choosing sequencing depth in a DNA storage system.

---

## References

R. Heckel, G. Mikutis, and R. N. Grass, "A Characterization of the DNA Data Storage Channel," *arXiv:1705.04732v1*, May 2017.