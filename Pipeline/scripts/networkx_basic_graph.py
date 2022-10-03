import itertools
import networkx as nx


def calculate_overlaps(p_a, p_b):
    ov = len(p_b)
    for i in range(min(len(p_a), len(p_b))):
        if p_a[-(i+1):] == p_b[:i+1]:
            ov = len(p_b) - (i+1)
    return ov


pep_count = int(snakemake.wildcards['pep_count'])
peptides = snakemake.params['pep'][:pep_count]

DG = nx.DiGraph()
DG.add_nodes_from(['Source', 'Sink'])

for peptide_a, peptide_b in itertools.product(peptides, peptides):
    if peptide_a == peptide_b:
        continue
    overlap = calculate_overlaps(peptide_a, peptide_b)
    DG.add_edge(peptide_a, peptide_b, weight=overlap)

DG.add_edges_from([('Source', p, {'weight': len(p)}) for p in peptides])
DG.add_edges_from([(p, 'Sink', {'weight': 0}) for p in peptides])

nx.write_gpickle(DG, snakemake.output[0])
