import re
import glob
import json
import pandas as pd
import altair as alt
from ProcessData import read_peptides


data = {'Peptide Count': [], 'Objective Value': [], 'Embedding': []}
files = ['/home/sara/Documents/VaccinesProject/ivp/Pipeline/17Aug_mhcI_with_mutations/ilp/lp_out/27030_vaccine_ilp_hog.json'] + \
        glob.glob('/home/sara/Documents/VaccinesProject/ivp/Pipeline/17Aug_mhcI_*embedded/ilp/lp_out/*_vaccine_ilp_hog.json')
for file in files:
    js = json.load(open(file))
    tmp = pd.DataFrame(js['SolutionInfo'])
    obj_val = tmp['ObjVal'][0]

    embed_len = re.search(r'17Aug_mhcI_(\d*)(embedded)*', file).group(1)
    if not embed_len:
        embed_len = 0
    else:
        embed_len = int(embed_len)

    pep_file = file.replace('lp_out', 'pep_out').replace('vaccine_ilp_hog.json', 'chosen_peptides_hog_inc_substrings.txt')
    peptides = read_peptides.main(pep_file)
    pep_count = len(peptides)

    data['Peptide Count'].append(pep_count)
    data['Objective Value'].append(obj_val)
    data['Embedding'].append(embed_len)

print(data)
df = pd.DataFrame(data)

sorting = ['mhcI', 'mhcI_1embedded', 'mhcI_2embedded', 'mhcI_3embedded', 'mhcI_4embedded', 'mhcI_5embedded', 'mhcI_6embedded', 'mhcI_7embedded', 'mhcI_8embedded', 'mhcI_9embedded', 'mhcI_10embedded']
domain = [df.min()['Objective Value'], (df.min()['Objective Value'] + df.max()['Objective Value'])/2, df.max()['Objective Value']]
print(domain)
# colors = ['#9AD6E5', '#66A6D9', '#BBB465', '#7A9B3A', '#DA0034', '#3A7A9B', '#FFAD00', '#9B3A7A', '#294B5D', '#DA2300', '#9B5B3A']
colors = ['#7A9B3A', '#DA0034', '#3A7A9B']
chart = alt.Chart(df).mark_circle(size=90).encode(
    x = alt.X('Embedding:Q'),
    y = alt.Y('Peptide Count'),
    # y = alt.Y('Objective Value', scale=alt.Scale(domain=[2, 2.8])),
    color = alt.Color('Objective Value', scale=alt.Scale(domain=domain, range=colors))
)

chart.save('cleavage_plots/embedding_vs_pepcount_vs_score_mhc1.html')