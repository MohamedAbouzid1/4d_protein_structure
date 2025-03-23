from functools import partial
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.features.nodes.amino_acid import amino_acid_one_hot

from ppi.edges import add_heavy_atom_connection
from ppi.graphs import prune_graph
from graphein.protein.graphs import construct_graph

def create_protein_graph(
    pdb_path: str,
    node_features: list = None,
    edge_threshold: float = 6.0,
    granularity: str = "atom",
    verbose: bool = False
):
    """
    Creates a protein graph from a PDB file with customizable features.
    
    Args:
        pdb_path (str): Path to PDB file
        node_features (list): List of node feature functions
        edge_threshold (float): Distance threshold for edge construction (Å)
        granularity (str): Graph granularity ('atom' or 'residue')
        verbose (bool): Whether to show verbose output
    
    Returns:
        ProteinGraph object
    """
    # Set default node features if not provided
    if node_features is None:
        node_features = [amino_acid_one_hot]
    
    # Create configuration
    config = ProteinGraphConfig(
        granularity=granularity,
        verbose=verbose,
        alt_locs="min_occupancy",
        node_metadata_functions=node_features,
        edge_construction_functions=[
            partial(add_heavy_atom_connection, threshold=edge_threshold),
        ],
        graph_metadata_functions=[
            partial(prune_graph, include="CA"),
        ]
    )
    
    # Construct and return graph
    return construct_graph(config=config, path=pdb_path)

