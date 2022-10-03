import vaccine_ilp_og
import vaccine_ilp_hog
import vaccine_ilp_naive
import gzip
import pickle
import networkx as nx
import pandas as pd
import gurobipy as gp
from gurobipy import GRB


def binarize_entries(df, th):
    df[df < th] = 0.0
    df[df >= th] = 1.0
    return df


path = snakemake.params['path']
k = int(snakemake.params['k'])
populations = snakemake.params['population']
peptides_count = snakemake.wildcards['pep_count']
peptides = snakemake.params['pep'][:int(peptides_count)]
min_hits = int(snakemake.params['min_hits'])

print('Read frequencies')
f_data = pd.read_pickle(snakemake.input['frequency_data'])
type = snakemake.params['type']

print('Read ba data')
ba_threshold = float(snakemake.params['binding_affinity_threshold'])
ba_matrix = pd.read_pickle(snakemake.input['binding_affinity_data'])
print('Binarize')
ba_matrix = binarize_entries(ba_matrix, ba_threshold)
print('Calc shared alleles')
full_known_alleles = gp.tuplelist(allele for allele in f_data.keys() if allele in ba_matrix.keys())
print(full_known_alleles)

# if int(snakemake.wildcards['pep_count']) <= 30:
#     coloring = True
# else:
#     coloring = False
coloring = False

if snakemake.params['approach'] == 'og':
    cost_graph = nx.read_gpickle(snakemake.input['overlap_graph'])

    vaccine_ilp_og.solve_vaccine_problem(k=k,
                                         alleles=full_known_alleles,
                                         freq_vector=f_data,
                                         B_matrix=ba_matrix,
                                         peptides=peptides,
                                         pep_count=peptides_count,
                                         graph=cost_graph,
                                         populations=populations,
                                         path=path,
                                         subtour_el='boeckler',
                                         coloring=coloring)

elif snakemake.params['approach'] == 'hog':
    print('Start HOGVAX')
    hog = nx.read_gpickle(snakemake.input['hog'])
    with open(snakemake.input['leaves'], 'rb') as handle:
        leaves_dict = pickle.load(handle)

    vaccine_ilp_hog.solve_msks_hog(k=k,
                                   alleles=full_known_alleles,
                                   freq_vector=f_data,
                                   B_matrix=ba_matrix,
                                   leaves=leaves_dict,
                                   pep_count=peptides_count,
                                   graph=hog,
                                   min_hits=min_hits,
                                   populations=populations,
                                   path=path,
                                   logging=False,
                                   coloring=coloring)

elif snakemake.params['approach'] == 'concat':
    print('Start concat ilp')
    vaccine_ilp_naive.solve_vaccine_problem(k=k,
                                            peptides=peptides,
                                            pep_count=peptides_count,
                                            alleles=full_known_alleles,
                                            freq_vector=f_data,
                                            B_matrix=ba_matrix,
                                            populations=populations,
                                            path=path)

else:
    print('Method is unclear, choose between og and hog for overlap graph or hierarchical overlap graph, respectively.')


print('Possible maximum for population %s: %g' % (populations, sum(f_data[i][pop] for pop in populations for i in f_data if i in ba_matrix.keys())))
# calculate maximum for given set of peptides --> if take all peptides from input set
# frequencies = set()
# for peptide in peptides:
#    for pop in populations:
#        for locus in full_known_alleles:
#            hit = f_data[locus][pop] * ba_matrix[locus][peptide]
#            if hit > 0:
#                frequencies.add((locus, f_data[locus][pop]))
#
# total_max_peptides = sum(f for _, f in frequencies)
# print('Possible maximum for given set of peptides: %g' % total_max_peptides)
