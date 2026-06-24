"""
Functions for manipulating river network data as Python NetworkX Graphs.
Matthew Heberger, April 2024

Functions to create a graph object, calculate the stream order (Strahler and Shreve)
# These are useful for visualizing and processing the river network data.

"""

import networkx as nx
from pandas import DataFrame


def make_river_network(df: DataFrame, terminal_node=None) -> nx.DiGraph:
    """
    Creates a network graph of our river network data.

    Arguments:
        df: a Pandas DataFrame with a unique id for each node in the network.
        The index should be a unique id for the node in the network. There
        must also there must be a field `nextdown`, the downstream neighbor,
        endcoding the direction of flow

        terminal_node: Optional: the id of the terminal node or most downstream node.

    Returns:
        A NetworkX DiGraph object.
    """
    G = nx.DiGraph()

    # Populate the graph's nodes and edges.
    for node_id, nextdown in df['nextdown'].items():
        # Add node with comid as node ID
        G.add_node(node_id)
        G.nodes[node_id]['area'] = df.at[node_id, 'unitarea']
        # Add edge from comid to nextdown
        if str(nextdown) != '0' and node_id != terminal_node:
            G.add_edge(node_id, nextdown)

    return G


def calculate_strahler_stream_order(graph: nx.DiGraph) -> nx.DiGraph:
    """
        Calculates the Strahler stream order for a river network that is
        represented as a NetworkX acyclic directed graph.
        Adds the attribute 'strahler_order' for nodes in the network.

        Arguments:
            a networkX acyclic directed Graph
        Returns:
            same Graph, with additional attribute on all nodes.

        """
    # Step 1: Determine upstream and downstream nodes
    up_nodes = {node: set() for node in graph.nodes()}
    down_nodes = {node: set() for node in graph.nodes()}
    for u, v in graph.edges():
        down_nodes[u].add(v)
        up_nodes[v].add(u)

    # Step 2: Assign initial stream orders
    for node in graph.nodes():
        graph.nodes[node]['strahler_order'] = 1

    # Step 3: Iterate through the nodes and update stream orders
    for node in nx.topological_sort(graph):
        upstream_orders = [graph.nodes[upstream]['strahler_order'] for upstream in up_nodes[node]]
        max_order = max(upstream_orders) if upstream_orders else 0
        order_count = upstream_orders.count(max_order)
        if order_count == 1:
            graph.nodes[node]['strahler_order'] = max_order
        else:
            graph.nodes[node]['strahler_order'] = max_order + 1

    return graph


def calculate_shreve_stream_order(graph: nx.DiGraph) -> nx.DiGraph:
    """
    Calculates the Shreve stream order for a river network that is
    represented as a NetworkX acyclic directed graph.
    Adds the attribute 'shreve_order' for nodes in the network.

    Arguments:
        a NetworkX DiGraph
    Returns:
        a NetworkX DiGraph with a new attribute 'shreve_order' for each node.

    """
    up_nodes = {node: set() for node in graph.nodes()}
    down_nodes = {node: set() for node in graph.nodes()}
    for u, v in graph.edges():
        down_nodes[u].add(v)
        up_nodes[v].add(u)

    for node in graph.nodes():
        graph.nodes[node]['shreve_order'] = 1

    for node in nx.topological_sort(graph):
        if up_nodes[node]:
            graph.nodes[node]['shreve_order'] = max([graph.nodes[upstream]['shreve_order']
                                                     for upstream in up_nodes[node]]) + 1

    return graph


def calculate_num_incoming(G: nx.DiGraph) -> nx.DiGraph:
    """calculates the number of upstream or incoming edges to each node
    and puts that in an attribute"""

    for node in G.nodes():
        num_incoming = G.in_degree(node)
        G.nodes[node]['num_incoming'] = num_incoming

    return G


def insert_node(G: nx.DiGraph, node, comid) -> nx.DiGraph:
    """
    Custom function to insert a new node in my flow network graph at a given location.

    To insert new river reaches and unit catchments into the existing flow network, it
    was helpful to create a Python NetworkX Graph.

    Adding an outlet point means inserting a new node to the graph.
    We selectively add and remove edges to maintain the connectivity of the graph.

    The logic is different if we are inserting a node
    into a unit catchment that is a "leaf" node, or one with a Strahler order 1
    or into a "stem" node, one with Strahler order > 1.
    """

    # Get the Strahler order of the node we're inserting into.
    order = G.nodes[comid]['strahler_order']

    if order == 1:
        # LEAF NODE
        # Add the new node to the network
        G.add_node(node)
        # Set some attributes for the node that will be useful later
        G.nodes[node]['new'] = True
        G.nodes[node]['type'] = 'leaf'
        # Add the network connection (edge) from the new node to the target unit catchment with comid = comid.
        G.add_edge(node, comid)

    else:
        # STEM NODE
        # Step 1: Add the new node
        G.add_node(node)
        G.nodes[node]['new'] = True
        G.nodes[node]['type'] = 'stem'

        # Step 2: Find incoming edges and remove them.
        predecessors = list(G.predecessors(comid))
        incoming_edges = list(G.in_edges(comid))
        for u, v in incoming_edges:
            G.remove_edge(u, v)

        # Step 3: Add an edge from the new node to the comid (node it is being inserted upstream of)
        G.add_edge(node, comid)

        # Step #4: Add new incoming edges to the new node
        for predecessor in predecessors:
            G.add_edge(predecessor, node)

    return G


def prune_node(G: nx.DiGraph, node) -> nx.DiGraph:
    """
    Prunes (removes) a node from a directed acyclic graph (DAG) and reconnects
     its upstream and downstream neighbors.

    Parameters:

        G (networkx.DiGraph): The directed acyclic graph.

        node (any hashable type): The node to be pruned from the graph.

    Returns:

        networkx.DiGraph: The DAG with the node pruned and neighbors reconnected.
    """
    if not G.has_node(node):
        raise ValueError("The specified node is not in the graph.")

    predecessors = list(G.predecessors(node))
    successors = list(G.successors(node))

    # Reconnect predecessors to successors
    for pred in predecessors:
        for succ in successors:
            if not G.has_edge(pred, succ):
                G.add_edge(pred, succ)

    # Remove the node from the graph
    G.remove_node(node)

    return G


def upstream_nodes(G: nx.DiGraph, node) -> list:
    """
    Return all upstream nodes for a given node in a directed acyclic graph

    Parameters:
        G (networkx.DiGraph): The directed acyclic graph.

        node (any hashable type): The node for which to find all upstream nodes.

    Returns:
        A list of upstream nodes.
    """
    if not G.has_node(node):
        raise ValueError("The specified node is not in the graph.")

    up_nodes = set()
    to_visit = [node]

    while to_visit:
        current = to_visit.pop()
        for predecessor in G.predecessors(current):
            if predecessor not in up_nodes:
                up_nodes.add(predecessor)
                to_visit.append(predecessor)

    return list(up_nodes)
