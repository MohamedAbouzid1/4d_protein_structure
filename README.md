# Protein Graph Analysis with Graphein and ORCA

This repository contains a set of tools to analyze protein structures using graph theory. It converts protein structures from PDB files into network graphs using Graphein and performs graphlet analysis using the ORCA tool.

## Overview

This project enables:
1. Converting protein structures to graph representations
2. Computing graphlet signatures using ORCA
3. Batch processing of multiple proteins from a FASTA file

## Installation

### Prerequisites

- Python 3.7 or higher
- NetworkX
- Graphein
- ORCA executable (must be compiled separately)

### Installing dependencies

```bash
pip install graphein networkx
```

For ORCA, you'll need to compile it from source:
- Download ORCA from: [https://github.com/thobye/orca](https://github.com/thobye/orca)
- Follow the compilation instructions in the ORCA repository

## Usage

### Command-line interface

```bash
python process_pdb_orca.py --pdb_dir /path/to/pdb/files --orca_path /path/to/orca/executable --fasta_file /path/to/your/sequences.fasta --output_dir /path/to/output
```

### Arguments

- `--pdb_dir`: Directory containing PDB files (expected format: "AF-{uniprot_id}-F1-model_v4.pdb")
- `--orca_path`: Path to the ORCA executable
- `--fasta_file`: Path to a FASTA file containing UniProt IDs
- `--output_dir`: Directory to store the ORCA results

## Components

### 1. Protein Graph Construction (`cust_construct_graph.py`)

This module uses Graphein to create protein graphs:

```python
from cust_construct_graph import create_protein_graph

# Create a protein graph from a PDB file
graph = create_protein_graph(
    pdb_path="path/to/protein.pdb",
    edge_threshold=6.0,  # in Angstroms
    granularity="atom"   # or "residue"
)
```

### 2. ORCA Interface (`prepare_and_run_orca.py`)

This module prepares graph data for ORCA and runs the graphlet analysis:

```python
from prepare_and_run_orca import prepare_and_run_orca

prepare_and_run_orca(
    g=graph,                         # NetworkX graph
    orca_path="/path/to/orca",       # Path to ORCA executable
    input_filename="input.in",       # ORCA input file
    output_filename="output.out",    # ORCA output file
    graphlet_size=5                  # Size of graphlets
)
```

### 3. Batch Processing (`process_pdb_orca.py`)

The main processing script that:
- Extracts UniProt IDs from a FASTA file
- Constructs graphs for each protein
- Runs ORCA analysis on each graph
- Saves results to the specified output directory

## ORCA Output Format

The ORCA tool produces files in the following format:
- Each line represents a node in the graph
- Each column represents a specific graphlet orbit count
- For graphlet size 5, there are 73 orbits

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- [Graphein](https://github.com/a-r-j/graphein) for protein graph creation
- [ORCA](https://github.com/thobye/orca) for graphlet counting

## Author

Mohamed Abouzid, 2025
