import sys
import networkx as nx
import read_peptides

# lp output
sol_file = snakemake.input['sol']
# og
og_file = snakemake.input['og']
# ilp chosen pep output
pep_file = snakemake.input['peptides']

edges = []
with open(sol_file) as f:
    for line in f:
        if line.startswith('#'):
            continue
        l = line.strip('\n').split(' ')
        # filter out hit variables named by alleles; DRB1 is not labelled with HLA prefix
        if 'HLA' in l[0] or 'DRB1' in l[0]:
            continue
        # filter for traversed edges
        if round(float(l[1])) == 0:
            continue
        edges.append(l[0])

edge_dict = {}
og = nx.read_gpickle(og_file)
for edge in og.edges():
    key = '_'.join(edge)
    edge_dict[key] = edge

# create subgraph from edges chosen by ILP in order to find eulerian path to construct vaccine sequence
sub_og = nx.DiGraph()
sub_og.add_edges_from([edge_dict[e] for e in edge_dict if e in edges])
print(sub_og)
print(nx.is_eulerian(sub_og))
seq = ''
if nx.has_eulerian_path(sub_og):
    print('Has Eulerian path')
    for e in nx.eulerian_path(sub_og):
        print(e)
        node_a, node_b = e
        length = og.edges[e]['weight']
        seq += node_b[len(node_b)-length:]
else:
    print('Attention! No Eulerian path!!')
    exit(-1)

print(seq)
print(len(seq))

chosen_peptides = read_peptides.main(pep_file)
for pep in chosen_peptides:
    if not pep in seq:
        print('{} does not occur in overall sequence!'.format(pep))
        exit(-1)

with open(snakemake.output[0], 'w') as out:
    out.write('> HOGVAX MHC optimized combined peptide vaccine sequence with overlaps\n')
    out.write(seq + '\n')
    out.write('> Concat MHC optimized combined peptide vaccine sequence concatenated\n')
    out.write(''.join(chosen_peptides))
