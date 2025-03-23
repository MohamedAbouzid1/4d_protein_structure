import networkx as nx
import subprocess

def prepare_and_run_orca(g, orca_path, input_filename, output_filename, graphlet_size=5):
    """
    Prepares the ORCA input file from a NetworkX graph, removes self-loops, and runs the ORCA tool.

    Parameters:
    - g1: NetworkX graph
    - orca_path: Path to the ORCA executable
    - input_filename: Name of the ORCA input file (default: 'orca_input.in')
    - output_filename: Name of the ORCA output file (default: 'orca_output.out')
    - graphlet_size: Size of the graphlets for ORCA (default: 5)
    """
    # Ensure g1 is a NetworkX graph
    if not isinstance(g, nx.Graph):
        raise TypeError("g1 is not a NetworkX graph.")

    # Remove self-loops
    self_loops = list(nx.selfloop_edges(g))
    g.remove_edges_from(self_loops)

    # Create a mapping from original node IDs to consecutive integers
    nodes = list(g.nodes())
    edges = list(g.edges())
    print(f"No. of nodes and edges after removing self nodes: {len(nodes)} nodes and {len(edges)} edges")
    node_mapping = {node: idx for idx, node in enumerate(nodes)}

    # Apply the mapping to edges
    mapped_edges = [(node_mapping[u], node_mapping[v]) for u, v in edges]

    # Write the ORCA input file
    with open(input_filename, "w") as f:
        f.write(f"{len(nodes)} {len(edges)}\n")
        for u, v in mapped_edges:
            f.write(f"{u} {v}\n")

    # Run the ORCA tool
    try:
        subprocess.run([orca_path, str(graphlet_size), input_filename, output_filename], check=True)
        print(f"ORCA analysis completed. Results saved in {output_filename}.")
    except subprocess.CalledProcessError as e:
        print(f"Error running ORCA: {e}")
