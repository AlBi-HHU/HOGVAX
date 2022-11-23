import networkx as nx


def draw_hog(hog, leaves, vals, key, draw_all=False):
    if draw_all:
        SG = hog
    else:
        SG = hog.edge_subgraph([(node_a, node_b) for node_a, node_b in vals if round(vals[node_a, node_b], 0) == 1.0])

    dot = nx.drawing.nx_pydot.to_pydot(SG)
    for i, e in enumerate(dot.get_edges()):
        attr = e.get_attributes()
        e.set_label(attr['weight'])
        if attr['is_slink'] == 'yes':
            e.set_color('blue')
            e.set_style('dashed')

        node_a, node_b = e.obj_dict['points']
        if round(vals[node_a, node_b], 0) == 1.0:
            e.set_color('red')
            if node_a in leaves.values():
                dot.get_node(node_a)[0].set_style('filled')
            if node_b in leaves.values():
                dot.get_node(node_b)[0].set_style('filled')
            dot.get_node(node_a)[0].set_color('red')
            dot.get_node(node_b)[0].set_color('red')

    for i, node in enumerate(dot.get_nodes()):
        attr = node.get_attributes()
        node.set_label(node.get_name() + ', ' + attr['string'] + ', ' + attr['interval_size'])

    dot.write_png(key)