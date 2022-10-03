import sys
import pickle
import pydot
import networkx as nx


def log(string):
    if logging_enabled:
        print(string)


def get_next():
    global counter
    counter += 1
    return str(counter)


def get_slink(c, node, graph):
    cur_node = node
    # loop until slink found and we return
    while True:
        for i in graph.successors(cur_node):
            if graph.nodes[i]['char'] == c:
                return i
        # if we are at root node and root has no child labeled with c, we return root node as slink
        if cur_node == '0':
            return cur_node
        # follow the suffix links further
        else:
            cur_node = graph.nodes[cur_node]['slink']


def add_slinks(graph):
    queue = ['0']
    while queue:
        cur = queue.pop(0)
        for i in graph.successors(cur):
            queue.append(i)
            # children of root node have suffix link back to root
            if cur == '0':
                graph.nodes[i]['slink'] = '0'
            # for all other nodes, follow suffix link of parent and check for child with same label as node i
            else:
                pred_slink = graph.nodes[cur]['slink']
                char = graph.nodes[i]['char']
                graph.nodes[i]['slink'] = get_slink(char, pred_slink, graph)


def get_or_add_node(c, prefix, node, graph):
    for i in graph.successors(node):
        if graph.nodes[i]['char'] == c:
            return i
    next_node = get_next()
    graph.add_node(next_node, char=c, string=prefix)
    graph.add_edge(node, next_node, weight=1)
    return next_node


def calc_lps(word):
    lps_table = [0] * len(word)
    lps_len = 0 # longest prefix suffix length is 0 at beginning
    i = 1   # first entry is always 0 so we start at index 1
    while i < len(word):
        if word[i] == word[lps_len]:
            lps_len += 1
            lps_table[i] = lps_len
            i += 1
        else:
            if lps_len > 0:
                lps_len = lps_table[lps_len-1]
            else:
                i += 1
    return lps_table


def insert(word, graph):
    cur_node = '0'
    word_nodes = [cur_node]
    lps = calc_lps(word)
    prefix = ''
    for c, l in zip(word, lps):
        prefix += c
        cur_node = get_or_add_node(c, prefix, cur_node, graph)
        word_nodes.append(cur_node)
        graph.nodes[cur_node]['ancestor'] = word_nodes[l]
    return cur_node


def main(count, strings, path, graph_out, leaves_out, logging=False, dot_bool=False):
    global logging_enabled
    logging_enabled = logging

    global counter
    counter = -1

    DG = nx.DiGraph()
    root = get_next()
    # we calculate ancestor here already for linear time hog construction
    DG.add_node(root, string='', ancestor='')

    log('Build trie')
    leaves = {}
    for i, entry in enumerate(strings):
        word = entry
        log('Insert word ' + word)
        leaves[word] = insert(word, DG)
    log('Built trie')

    log('Add slinks')
    add_slinks(DG)

    nx.write_gpickle(DG, graph_out)
    with open(leaves_out, 'wb') as handle:
        pickle.dump(leaves, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # draw graph
    if dot_bool:
        print('Write dot')
        dot = nx.drawing.nx_pydot.to_pydot(DG)

        for i, node in enumerate(dot.get_nodes()):
            attr = node.get_attributes()
            # use int here cause root is first node in list and it is easier to check for index here
            if i == 0:
                node.set_label('0 root')
            else:
                node.set_label(node.get_name() + ', ' + attr['string'])
                edge = pydot.Edge(node.get_name(), attr['slink'], color='blue', style='dashed')
                dot.add_edge(edge, )

        dot.write_png(path + count + '_aho_corasick.png')


if __name__ == '__main__':
    count = int(snakemake.wildcards['pep_count'])
    strings = snakemake.params['pep'][:count]
    path = snakemake.params['path']
    logging = snakemake.params['log']
    out1 = snakemake.output[0]
    out2 = snakemake.output[1]

    # if count <= 30:
    #     coloring = True
    # else:
    #     coloring = False
    # TODO
    coloring = False

    main(str(count), strings, path, out1, out2, logging=logging, dot_bool=coloring)
