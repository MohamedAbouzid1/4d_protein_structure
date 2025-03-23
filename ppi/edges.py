from graphein.protein.edges.distance import (
    compute_distmat, get_interacting_atoms
)
from graphein.protein.utils import filter_dataframe
import networkx as nx

from ppi.utils import get_ca_atom


def add_heavy_atom_connection(
    G: nx.Graph,
    threshold: float,
) -> nx.Graph:
    """Connects the CA atoms of residues if any heavy atom lies within the
    threshold.

    Parameters
    ----------
    G : nx.Graph
        Graph to add edges to
    threshold : float
        Threshold distance for heavy atom interaction in Ångstrom

    Returns
    -------
    nx.Graph
        Graph with added edges
    """
    if 'atomic_dist_mat' in G.graph:
        dist_mat = G.graph['atomic_dist_mat']
    else:
        pdb_df = filter_dataframe(
            dataframe=G.graph["pdb_df"],
            by_column="node_id",
            list_of_values=list(G.nodes()),
            boolean=True,
        )
        dist_mat = compute_distmat(pdb_df)

    interacting_nodes = get_interacting_atoms(threshold, dist_mat)
    interacting_nodes = zip(interacting_nodes[0], interacting_nodes[1])

    for a1, a2 in interacting_nodes:
        df = G.graph["pdb_df"]
        node_1_name = df.loc[df.index[a1], "node_id"]
        node_2_name = df.loc[df.index[a2], "node_id"]

        ca_1 = get_ca_atom(G, node_1_name)
        ca_2 = get_ca_atom(G, node_2_name)

        if ca_1 == '' or ca_2 == '':
            continue

        # We are in the same residue, so we skip
        if ca_1 == ca_2:
            continue

        # If a heavy atom was already within the threshold, there is no need
        # to connect the CA atoms again.
        if G.has_edge(ca_1, ca_2):
            continue

        if ca_1 not in G.nodes or ca_2 not in G.nodes:
            raise RuntimeError(
                "An error occured while trying to obtain the correct CA atoms."
            )

        G.add_edge(ca_1, ca_2, kind={"heavy_atom_connection"})

    return G
