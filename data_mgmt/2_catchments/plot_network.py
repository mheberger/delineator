import networkx as nx
import graphviz

SHOW_AREA = True
VERTICAL_PLOT = True
DIAGRAM_FORMAT = 'png'


def area_to_size(area, max_area):
    if max_area is None:
        return 10
    else:
        size = 0.2 + 1.2 * (area/max_area) ** 0.5
    return size


def draw_graph(G: nx.DiGraph, filename: str, title="River Network Graph"):
    """
    Plots a NetworkX directed acyclic graph (dag)
    using the GraphViz library. Note that you need to install this software and it needs to be
    on the system path.

    Args:
        G: the river network directed graph.
        title: Descriptive text for the title
        filename: the filename (bare, no extension) to create
        vertical: if True, will arrange plot from top to bottom,
          Default is left to right.
    """

    # Since we are going to size the nodes by area, find the largest so we can scale appropriately
    max_area_node = max(G.nodes, key=lambda node: G.nodes[node].get('area', 0))
    try:
        max_area = G.nodes[max_area_node]['area']
    except KeyError:
        max_area = None

    # Initialize a new directed graph in Graphviz
    if VERTICAL_PLOT:
        dot = graphviz.Digraph()
    else:
        dot = graphviz.Digraph(graph_attr={'rankdir': 'LR'})

    dot.graph_attr['label'] = title
    dot.graph_attr['labelloc'] = 't'  # 't' for top, 'b' for bottom, 'c' for center
    dot.graph_attr['fontsize'] = '20'  # Set the font size for the title

    # Add nodes with distance attributes as labels
    for node in G.nodes():
        label = f"{node}"

        if SHOW_AREA:
            if 'area' in G.nodes[node]:
                area = round(G.nodes[node]['area'])
                label = f"{label}\n{area}"
                size = area_to_size(G.nodes[node]['area'], max_area)
                size = str(round(size, 1))
                dot.node(str(node), width=size, height=size, fixedsize='true')

        if G.nodes[node].get('custom'):
            dot.node(str(node), label=label, style='filled', fillcolor='lightblue', fontname='Arial')
        else:
            dot.node(str(node), label=label, fontname='Arial')

    # Add edges
    for u, v in G.edges():
        dot.edge(str(u), str(v))

    # Save the graph as a file and render it
    dot.render(filename, format=DIAGRAM_FORMAT, cleanup=True)

    # Display the graph using PIL
    # image = Image.open(filename + "." + DIAGRAM_FORMAT)
    # image.show()


if __name__ == "__main__":
    # Example usage:
    # Create a directed acyclic graph
    GR = nx.DiGraph()
    GR.add_edges_from([(1, 2), (1, 3), (3, 4), (2, 4), (4, 5)])

    # Plot the DAG
    draw_graph(GR, "my_dag")
