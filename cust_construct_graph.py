import graphein.protein as gp
from functools import partial
from graphein.protein.features.nodes.amino_acid import amino_acid_one_hot
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
        node_features = [gp.amino_acid_one_hot]
    
    # Define edge construction functions
    dist_edge_func = {
        "edge_construction_functions": [
            partial(gp.add_distance_threshold, threshold=edge_threshold, long_interaction_threshold=0)
        ]
    }
    
    # Define node metadata functions
    one_hot = {
        "node_metadata_functions": node_features
    }
    
    # Create configuration by combining components
    config = gp.ProteinGraphConfig(
        granularity=granularity,
        verbose=verbose,
        alt_locs="min_occupancy",
        **{**dist_edge_func, **one_hot}
    )
    
    # Construct and return graph
    return construct_graph(config=config, path=pdb_path)