import os
import glob
import json
import re
import gzip
import pandas as pd
import altair as alt
from EvaluateVaccines import evaluate_vaccine as ev
from EvaluateVaccines import comparison_data as cd
from ProcessData import read_peptides


mhc1_all_ba = pd.read_pickle(gzip.open('../../Gifford_Data/25June_mhc1_netmhc-4.1_pred_affinity_pivot.pkl.gz'))
mhc1_af = pd.read_pickle('../../Gifford_Data/IEDB_population_frequency2392_normalized.pkl')

mhc2_all_ba = pd.read_pickle(gzip.open('../../Gifford_Data/25June_mhc2_netmhcii-4.1_pred_affinity_pivot_v1v2.pkl.gz'))
mhc2_af = pd.read_pickle('../../Gifford_Data/IEDB_population_frequency_mhc2_275normalized.pkl')

mhc1_scores = {'Size': [], 'Score': [], 'Method': []}
for file in glob.glob('/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/optivax_1-19_vaccines_mhc1/13Sept_optivax_unlinked_peptides_*.txt'):
    size = re.search(r'13Sept_optivax_unlinked_peptides_(\d+)', file).group(1)
    mhc1_scores['Size'].append(size)

    peptides = read_peptides.main(file)
    mhc1_scores['Score'].append(cd.get_cmp_obj(peptides, mhc1_all_ba, mhc1_af, smaller_to_zero=True).pred_at_least_one('World'))
    mhc1_scores['Method'].append('OptiVax')

print(mhc1_scores)

loci_cov = {'Locus': [], 'Score': [], 'Size': []}
for hog_file in glob.glob('/home/sara/Documents/VaccinesProject/ivp/Pipeline/15Sept_mhcI_*/'):
    size = re.search(r'15Sept_mhcI_(\d+)', hog_file).group(1)
    mhc1_scores['Size'].append(size)

    peptides = read_peptides.main(hog_file + 'ilp/pep_out/4512_chosen_peptides_hog_inc_substrings.txt')
    mhc1_scores['Score'].append(cd.get_cmp_obj(peptides, mhc1_all_ba, mhc1_af, smaller_to_zero=True).pred_at_least_one('World'))
    mhc1_scores['Method'].append('HOGVAX')

print(mhc1_scores)

df = pd.DataFrame.from_dict(mhc1_scores)
print(df)

chart = alt.Chart(df).mark_line().encode(
    x='Size:Q',
    y=alt.Y('Score:Q', scale=alt.Scale(domain=[0.4, 1.0])),
    color=alt.Color('Method:N', scale=alt.Scale(range=['#DA0034', '#FFC12C']))
)
chart.save('plots/evalvax_cmp_mhc1.html')


mhc2_scores = {'Size': [], 'Score': [], 'Method': []}
for file in glob.glob('/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/optivax_1-19_vaccines_mhc2/13Sept_optivax_unlinked_peptides_*.txt'):
    size = re.search(r'13Sept_optivax_unlinked_peptides_(\d+)', file).group(1)
    mhc2_scores['Size'].append(size)

    peptides = read_peptides.main(file)
    mhc2_scores['Score'].append(cd.get_cmp_obj(peptides, mhc2_all_ba, mhc2_af, smaller_to_zero=True).pred_at_least_one('World'))
    mhc2_scores['Method'].append('OptiVax')

print(mhc2_scores)

for hog_file in glob.glob('/home/sara/Documents/VaccinesProject/ivp/Pipeline/22Sept_mhc2_*/'):
    size = re.search(r'22Sept_mhc2_(\d+)', hog_file).group(1)
    mhc2_scores['Size'].append(size)

    pepfile = hog_file + 'ilp/pep_out/37435_chosen_peptides_hog_incl_substrings.txt'
    if not os.path.isfile(pepfile):
        pepfile = hog_file + 'ilp/pep_out/37435_chosen_peptides_hogvax_incl_substrings.txt'
    peptides = read_peptides.main(pepfile)
    mhc2_scores['Score'].append(cd.get_cmp_obj(peptides, mhc2_all_ba, mhc2_af, smaller_to_zero=True).pred_at_least_one('World'))
    mhc2_scores['Method'].append('HOGVAX')

print(mhc2_scores)

df = pd.DataFrame.from_dict(mhc2_scores)
print(df)

chart = alt.Chart(df).mark_line().encode(
    x='Size:Q',
    y=alt.Y('Score:Q', scale=alt.Scale(domain=[0.6, 1.0])),
    color=alt.Color('Method:N', scale=alt.Scale(range=['#DA0034', '#FFC12C']))
)
chart.save('plots/evalvax_cmp_mhc2.html')