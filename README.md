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
