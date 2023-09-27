from collections import defaultdict

def aho_corasick_algorithm(tree, leaves, peptides):
    # sort peptides in reversed order
    sorted_peptides = sorted(peptides, key=len, reverse=True)
    superstrings = defaultdict(list)
    for i, peptide in enumerate(sorted_peptides):
        if i % 100 == 0:
            print(i, 'of', len(sorted_peptides))
        # stop when reaching peptides of the shortest length
        if len(peptide) == len(sorted_peptides[-1]):
            break
        # thread characters through trie starting at root node
        state = '0'
        for c in peptide:
            while c not in [tree.nodes[succ]['char'] for succ in tree.successors(state)]:
                if state == '0':
                    return False
                else:
                    state = tree.nodes[state]['slink']
            # get next state
            state = [s for s in tree.successors(state) if tree.nodes[s]['char'] == c][0]
            if tree.nodes[state]['string'] in leaves:
                superstrings[peptide].append(tree.nodes[state]['string'])
            if 'dictlink' in tree.nodes[state]:
                superstrings[peptide] = superstrings[peptide] + tree.nodes[state]['dictlink']

    return superstrings
