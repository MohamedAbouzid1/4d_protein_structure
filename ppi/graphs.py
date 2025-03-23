import networkx as nx

from ppi.utils import get_ca_atom


def heavy_atom_coords_to_ca(G: nx.Graph) -> nx.Graph:
    """Adds all heavy atom coordinates to the alpha carbon atom of the residue.
    This is required to later calculate the labels between to graphs.
    """
    heavy_atom_key = 'heavy_atom_coords'
    for node_name in G.nodes:
        ca_atom_name = get_ca_atom(G, node_name)

        if ca_atom_name == '':
            continue

        if node_name == ca_atom_name:
            continue

        ca_atom = G.nodes[ca_atom_name]
        atom = G.nodes[node_name]

        # Only add the heavy atom coordinates if the atoms are from the same
        # residue.
        if ca_atom['residue_name'] != atom['residue_name']:
            continue

        if heavy_atom_key not in ca_atom:
            ca_atom[heavy_atom_key] = []
            ca_atom[heavy_atom_key].extend(ca_atom['coords'])

        ca_atom[heavy_atom_key].extend(G.nodes[node_name]['coords'])

    return G


def prune_graph(
    G: nx.Graph,
    include: str | list[str] = 'CA',
) -> nx.Graph:
    """Prunes the graph to only include the specified atoms.

    Parameters
    ----------
    G : nx.Graph
        Graph to prune
    include : str | list[str], optional
        Atom types that should be included in the graph, by default 'CA'

    Returns
    -------
    nx.Graph
        Pruned graph
    """

    if include == 'CA':
        include = ['CA']

    nodes_to_remove = [
        node for node in G.nodes() if node.split(":")[-1] not in include
    ]

    G.remove_nodes_from(nodes_to_remove)

    return G
