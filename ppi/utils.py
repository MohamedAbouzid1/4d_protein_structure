__author__ = ['Michael Variola']
__email__ = ['michael.variola@gmail.com']
import functools
import pandas as pd
import networkx as nx

@functools.lru_cache(maxsize=None)
def get_atom_df(
    G: nx.Graph,
    chain_id,
    residue_number,
    atom_type,
) -> pd.DataFrame:
    df = G.graph['pdb_df']
    return df.loc[
        (df['chain_id'] == chain_id) &
        (df['residue_number'] == residue_number) &
        (df['atom_name'] == atom_type)
    ]


def get_ca_atom(G: nx.Graph, node_name: str) -> str:

    # First we try a fast search based on the name of the node, before we start
    # slow lookups in the dataframe
    try:
        # The last part _ is the atom_type, which is not required.
        chain, residue_name, residue_number, alt_loc, _ = node_name.split(':')
        alt_loc = f":{alt_loc}"
    except ValueError:
        # At this point the node has no alt_loc, the other parts are always
        # available
        alt_loc = ''
        chain, residue_name, residue_number, _ = node_name.split(':')

    possible_ca_name = f"{chain}:{residue_name}:{residue_number}{alt_loc}:CA"

    if possible_ca_name in G.nodes:
        return possible_ca_name

    # If the CA atom of the given atom does not match the pattern above, we
    # do a more thorough search which takes longer.
    node = G.nodes[node_name]
    # This is pretty fast and can easily be cached
    result = get_atom_df(G, node['chain_id'], node['residue_number'], 'CA')

    # In larger proteins there is a chance for ambiguety, if that is the case
    # we search further
    if len(result) > 1:
        result = result.loc[result['residue_name'] == node['residue_name']]
        # There is still the possibility to run into errors when insertions
        # and alternative-locations occur
        if len(result) > 1:
            node_df = G.graph['pdb_df']

            node_df = node_df.loc[node_df['node_id'] == node_name]

            if len(node_df) > 1:
                node_df = node_df.loc[
                    (node_df['x_coord'] == node['coords'][0]) &
                    (node_df['y_coord'] == node['coords'][1]) &
                    (node_df['z_coord'] == node['coords'][2])
                ]

            result = result.loc[
                (result['insertion'] == node_df['insertion'].item()) &
                (result['alt_loc'] == node_df['alt_loc'].item())
            ]

    if len(result) != 1:
        # At this point twe cannot do anything else. This can happen if an
        # atom has no corresponding CA atom
        return ''

    return result['node_id'].item()