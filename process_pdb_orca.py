import os
import argparse
from cust_construct_graph import create_protein_graph
from prepare_and_run_orca import prepare_and_run_orca

def process_pdb(uniprot_id, pdb_dir, orca_path, output_dir):
    """Processes a single PDB file with Graphein and ORCA."""
    pdb_filename = f"AF-{uniprot_id}-F1-model_v4.pdb"
    pdb_path = os.path.join(pdb_dir, pdb_filename)

    if not os.path.exists(pdb_path):
        print(f"PDB file not found: {pdb_path}")
        return

    try:
        # Construct graph
        g = create_protein_graph(pdb_path=pdb_path)
        print(f"Graph constructed for {uniprot_id} with {len(g.nodes())} nodes and {len(g.edges())} edges")

        # Prepare ORCA I/O paths
        os.makedirs(output_dir, exist_ok=True)
        input_file = os.path.join(output_dir, f"{uniprot_id}.in")
        output_file = os.path.join(output_dir, f"{uniprot_id}.out")

        # Run ORCA
        prepare_and_run_orca(g, orca_path, input_file, output_file)

    except Exception as e:
        print(f"Error processing {uniprot_id}: {e}")

def process_fasta(fasta_path, pdb_dir, orca_path, output_dir):
    """Processes all entries in the FASTA file."""
    if not os.path.exists(fasta_path):
        print(f"FASTA file not found: {fasta_path}")
        return
    
    # Extract UniProt IDs from FASTA
    uniprot_ids = []
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                uniprot_id = line[1:].split()[0]
                uniprot_ids.append(uniprot_id)

    # Process each UniProt ID
    for uid in uniprot_ids:
        print(f"Processing {uid} from {fasta_path}...")
        process_pdb(uid, pdb_dir, orca_path, output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process PDB structures and run ORCA.")
    parser.add_argument("--pdb_dir", type=str, required=True, help="Directory containing PDB files")
    parser.add_argument("--orca_path", type=str, required=True, help="Path to ORCA executable")
    parser.add_argument("--fasta_file", type=str, required=True, help="Path to the FASTA file")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to store ORCA results")

    args = parser.parse_args()

    print(f"Processing FASTA file: {args.fasta_file}")
    process_fasta(args.fasta_file, args.pdb_dir, args.orca_path, args.output_dir)
