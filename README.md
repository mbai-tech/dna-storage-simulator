# DNA Storage Simulator

A Python-based simulation system for storing information in synthetic DNA molecules.

## Overview

This simulator demonstrates how digital information can be encoded into DNA sequences using M molecules of length L nucleotides. Each nucleotide (A, T, G, C) represents 2 bits of information using quaternary encoding.

## Features

- **Flexible Configuration**: Specify any number of molecules (M) and molecule length (L)
- **Text Encoding/Decoding**: Convert text to DNA sequences and back
- **Storage Capacity Calculation**: Automatic calculation of storage capacity in bits, bytes, and KB
- **Error Simulation**: Simulate sequencing errors to test data integrity
- **Capacity Management**: Automatic validation against storage limits

## Storage Encoding

- Each nucleotide encodes 2 bits:
  - `A` = `00`
  - `T` = `01`
  - `G` = `10`
  - `C` = `11`

- Storage capacity: **M × L × 2 bits**
  - Example: 50 molecules × 100 nucleotides = 10,000 bits = 1,250 bytes

## Usage

### Basic Example

```python
from dna_simulator import DNAStorageSimulator

# Initialize with M=50 molecules of L=100 nucleotides each
simulator = DNAStorageSimulator(num_molecules=50, molecule_length=100)

# Store data
message = "Hello, DNA storage!"
molecules = simulator.store_data(message)

# Retrieve data
retrieved = simulator.retrieve_data()
print(retrieved)  # "Hello, DNA storage!"
```

### Get Storage Information

```python
info = simulator.get_storage_info()
print(f"Capacity: {info['capacity_bytes']} bytes")
print(f"Molecules used: {info['molecules_used']}")
```

### Simulate Errors

```python
# Simulate 5% error rate (typical for DNA sequencing)
noisy_molecules = simulator.simulate_errors(error_rate=0.05)
retrieved_noisy = simulator.retrieve_data(noisy_molecules)
```

### Display Molecules

```python
# Show first 10 molecules
simulator.display_molecules(max_display=10)
```

## Running the Demo

```bash
python demo.py
```

The demo includes:
1. Basic storage and retrieval
2. Error simulation
3. Different data sizes
4. Capacity limit testing

## File Structure

- `dna_simulator.py` - Main simulator class and encoding/decoding functions
- `demo.py` - Demonstration script with examples
- `README.md` - This file

---

## Recovery Simulator (`src/dna_storage_sim.py`)

A separate simulation tool that models the sequencing retrieval process and verifies recovery rates against theoretical predictions from [Heckel et al. 2017 (arXiv:1705.04732)](https://arxiv.org/abs/1705.04732).

### Setup

```bash
cd "dna-storage-simulator"
python3 -m venv venv
source venv/bin/activate
pip install numpy matplotlib
```

### Running

```bash
# Activate the virtual environment first
source venv/bin/activate

# Basic run (default 30 M values)
python src/dna_storage_sim.py

# With beta=2 and custom M range
python src/dna_storage_sim.py --beta 2.0 --M_values 1 2 5 10 100 1000 10000 100000

# With reliability check
python src/dna_storage_sim.py --beta 2.0 --target_recovery 0.60

# Find optimal coverage c* for a given cost ratio
python src/dna_storage_sim.py --beta 2.0 --target_recovery 0.60 --find_c_star --q 10000

# Full q-sweep with design recommendation plots
python src/dna_storage_sim.py --beta 2.0 --target_recovery 0.60 \
  --q_values 1000 10000 100000 --q_output q_sweep.csv

# Save simulation results to CSV
python src/dna_storage_sim.py --beta 2.0 --output results.csv
```

### CLI Options

| Flag | Default | Description |
|---|---|---|
| `--M_values` | 30 log-spaced points | List of M values (number of unique strands) |
| `--c` | 1.0 | Coverage depth (reads per unique strand) |
| `--beta` | 1.0 | Strand-length factor: L = beta × log₂(M); use > 1 |
| `--trials` | 50 | Independent simulation trials per M value |
| `--mode` | `per_unique` | `per_unique` (N=c×M) or `per_physical` (N=c×M×DUP) |
| `--q` | 10000 | Synthesis-to-sequencing cost ratio |
| `--target_recovery` | None | Reliability threshold: mean − 2×SE ≥ target |
| `--find_c_star` | off | Search for the cost-minimising coverage c* |
| `--q_values` | None | Sweep multiple q values and generate design plots |
| `--q_output` | `q_sweep.csv` | Output CSV for q-sweep results |
| `--output` | None | Output CSV for main simulation sweep |

### Model

- **M** unique strands, each duplicated **DUP=10** times before sequencing
- Each read samples a strand uniformly at random (equivalent to sampling a unique strand ID in [0, M))
- Fraction recovered: `Q/M` where Q = number of distinct strand IDs seen
- Theoretical limit (Heckel et al., Theorem 1): `1 − e^{−c}`

### Cost Model (Heckel et al. 2017)

| Symbol | Formula | Meaning |
|---|---|---|
| `Rs` | `(1 − e^{−c})(1 − 1/β)` | Bits stored per nucleotide synthesized |
| `Rr` | `Rs / c` | Bits recovered per nucleotide sequenced |
| `cost` | `q/Rs + 1/Rr` | Cost per stored bit |

For q = 10,000 (typical), the optimal coverage is c* ≈ 9.2.

### File Structure

```
dna-storage-simulator/
├── src/
│   ├── dna_storage_sim.py   # Simulation engine and CLI
│   └── cost_model.py        # Theoretical cost model (Rs, Rr, cost)
├── plots/                   # Output plots (PNG)
├── venv/                    # Python virtual environment
└── README.md
```

## Parameters

### DNAStorageSimulator(num_molecules, molecule_length)

- **num_molecules** (M): Number of DNA molecules available for storage
- **molecule_length** (L): Length of each molecule in nucleotides

### Methods

- `store_data(data: str)` - Store text data and return DNA molecules
- `retrieve_data(molecules: List[str])` - Decode molecules back to text
- `simulate_errors(error_rate: float)` - Add random errors to molecules
- `get_storage_info()` - Get capacity and usage statistics
- `display_molecules(max_display: int)` - Print DNA sequences

## Example Output

```
Original message: "Hello, DNA storage!"
Data stored in 1 molecules

Molecule 1: TGGACAATGCTATGGAGTGTTGATTGGATTATGGATATCCCTGGCGATTGATGGACATGA...

Retrieved message: "Hello, DNA storage!"
✓ SUCCESS: Retrieved message matches original!
```

## Notes

- Text is converted to binary (8 bits per character)
- Binary data is encoded as DNA using 2 bits per nucleotide
- Padding is automatically added to fill molecules
- Null characters are stripped during retrieval
- Capacity limits are enforced during storage
