import sys
import networkx as nx
import read_peptides

# lp output
sol_file = snakemake.input['sol']
# hog
hog_file = snakemake.input['hog']
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
hog = nx.read_gpickle(hog_file)
for edge in hog.edges():
    node_a, node_b = edge
    key = hog.nodes[node_a]['string'] + '_' + hog.nodes[node_b]['string']
    edge_dict[key] = edge

# create subgraph from edges chosen by ILP in order to find eulerian path to construct vaccine sequence
sub_hog = nx.DiGraph()
sub_hog.add_edges_from([edge_dict[e] for e in edge_dict if e in edges])
print(sub_hog)
print(nx.is_eulerian(sub_hog))
seq = ''
for e in nx.eulerian_circuit(sub_hog):
    if hog.edges[e]['is_slink'] == 'yes':
        continue
    else:
        node_a, node_b = e
        length = hog.edges[e]['weight']
        node_seq = hog.nodes[node_b]['string']
        seq += node_seq[len(node_seq)-length:]

print(seq)
print(len(seq))

chosen_peptides = read_peptides.main(pep_file)
for pep in chosen_peptides:
    if not pep in seq:
        print('{} does not occur in overall sequence!'.format(pep))
        exit(-1)

with open(snakemake.output[0], 'w') as out:
    out.write('> HOGVAX optimized combined peptide vaccine sequence with overlaps\n')
    out.write(seq + '\n')
    out.write('> Concat optimized combined peptide vaccine sequence concatenated\n')
    out.write(''.join(chosen_peptides))
