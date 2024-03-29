import os
import pydot
import pickle
import networkx as nx
import aho_corasick_trie
from collections import defaultdict


def log(string):
    if logging_enabled:
        print(string)


def contract_ac(unmarked_nodes, slinks, graph):
    for node in unmarked_nodes:
        parent = list(graph.pred[node])[0]
        for succ in graph.successors(node):
            w = len(graph.nodes[succ]['string']) - len(graph.nodes[parent]['string'])
            graph.add_edge(parent, succ, weight=w)
            graph.nodes[succ]['char'] = graph.nodes[node]['char']

        if node in slinks:
            for link_node in slinks[node]:
                if link_node in graph.nodes():
                    graph.nodes[link_node]['slink'] = graph.nodes[node]['slink']
        graph.remove_node(node)
    return graph


def get_incoming_slinks(unmarked_nodes, graph):
    incoming_slinks = {}
    for node in graph.nodes():
        # skip root node
        if node == '0':
            continue
        # slink is incoming edge to slink_suc
        slink_suc = graph.nodes[node]['slink']
        # skip nodes that will not be removed
        if slink_suc not in unmarked_nodes:
            continue

        if slink_suc in incoming_slinks:
            incoming_slinks[slink_suc].append(node)
        else:
            incoming_slinks[slink_suc] = [node]
    return incoming_slinks


def mark_nodes(leaves, graph):
    # root node is initially marked
    marked_nodes = {'0'}
    # for each node initialize child nodes as empty set
    children = defaultdict(lambda: set())
    for leaf in leaves.values():
        cur_node = leaf
        marked_nodes.add(cur_node)
        # while current node is not root node
        while cur_node != '0':
            interval = graph.nodes[cur_node]['interval_size']
            for child in children[cur_node]:
                interval -= graph.nodes[child]['interval_size']
            if interval > 0:
                marked_nodes.add(cur_node)
            children[cur_node] = set()
            children[graph.nodes[cur_node]['ancestor']].add(cur_node)
            cur_node = graph.nodes[cur_node]['slink']

    return marked_nodes


def compute_interval_size(node, graph):
    if graph.out_degree(node) == 0:
        graph.nodes[node]['interval_size'] = 0
        return 1
    else:
        num_leaves = sum(compute_interval_size(suc, graph) for suc in graph.successors(node))
        graph.nodes[node]['interval_size'] = num_leaves
        return num_leaves


def compute_hog(pep_count, leaves_dict, aho_corasick, path, logging=False, dot_bool=False):
    global logging_enabled
    logging_enabled = logging

    # compute |R(u)| for each node u starting with root = 0
    log('Compute interval sizes')
    compute_interval_size('0', aho_corasick)
    # mark nodes
    log('Execute linear time algorithm')
    md_nodes = mark_nodes(leaves_dict, aho_corasick)

    # get nodes that we have to remove from aho corasick to get hog
    rm_nodes = set(aho_corasick.nodes) - md_nodes
    log('Nodes to be removed' + ', '.join(rm_nodes))

    # need to temporary store suffix link nodes for correct contraction
    incoming_slinks = get_incoming_slinks(rm_nodes, aho_corasick)

    # remove nodes and contract edges
    log('Contract aho corasick')
    hog = contract_ac(rm_nodes, incoming_slinks, aho_corasick)

    # set suffix links as explicit edges -> need this for ILP subtour elimination constraints
    nx.set_edge_attributes(hog, 'no', name='is_slink')
    slinks = [(i, hog.nodes[i]['slink'], {'weight': 0, 'is_slink': 'yes'}) for i in hog.nodes if i != '0']
    hog.add_edges_from(slinks)

    # write hog
    if not os.path.exists(path + 'pickle_files/'):
        os.mkdir(path + 'pickle_files/')
    with open(path + 'pickle_files/' + pep_count + '_hog.gpickle', 'wb') as file:
        pickle.dump(hog, file, protocol=pickle.HIGHEST_PROTOCOL)

    if dot_bool:
        hog_copy = hog.copy()
        for edge in hog_copy.edges:
            hog_copy.edges[edge]['weight'] = str(hog_copy.edges[edge]['weight'])

        dot = nx.drawing.nx_pydot.to_pydot(hog_copy)
        for i, e in enumerate(dot.get_edges()):
            attr = e.get_attributes()
            e.set_label(attr['weight'])
            if attr['is_slink'] == 'yes':
                e.set_color('blue')
                e.set_style('dashed')

        for i, node in enumerate(dot.get_nodes()):
            attr = node.get_attributes()
            node.set_label(node.get_name() + ', ' + attr['string'] + ', ' + attr['interval_size'])

        if not os.path.exists(path + 'vis/'):
            os.mkdir(path + 'vis/')
        dot.write_png(path + 'vis/' + pep_count + '_hog.png')

    return leaves_dict, hog
