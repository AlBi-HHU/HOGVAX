"""
Compare cleavage prediction peptides of different (embedded) vaccines with cleavage predicted peptides from SARS-CoV-2
proteome.
"""
import re
import glob
import pandas as pd
import altair as alt
from EvaluateVaccines import calculate_cleaved_peptides


def read_peptides(file):
    peptides = []
    with open(file, 'r') as f:
        for line in f:
            pep = line.strip()
            peptides.append(pep)
    print(len(peptides))
    return peptides


def transform_data(data):
    new_dict = {'embedding': [], 'tool': [], 'fraction': []}
    for key in data:
        new_dict['embedding'].extend([key, key])
        new_dict['tool'].extend(['NetChop3.1', 'Pepsickle'])
        new_dict['fraction'].extend([data[key]['NetChop3.1'], data[key]['Pepsickle']])
    return new_dict


out_name = '23Aug_mhcI_cleaved_virus_proteins'
sars_peptides_pepsickle = read_peptides('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/SarsCoVCleavage/pepsickle_sars_cov_2_cleavage_peptides.txt')
sars_peptides_netchop = read_peptides('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/SarsCoVCleavage/netchop_sars_cov_2_cleavage_peptides.txt')

intersection_count = {}
files = glob.glob('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/23AugEmbeddingPredictions/MHCI/NetChop3.1_17Aug_mhcI_*_reformat.txt') + \
        ['/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/23AugEmbeddingPredictions/MHCI/NetChop3.1_17Aug_mhcI_with_mutations_reformat_concat.txt']
for file in files:
    match = re.search(r'_17Aug_(mhcI+(_\d+embedded)*(_hogvax)*)(_with_mutations)*(_reformat_concat)*', file)
    name = match.group(1).replace('_hogvax', '')
    if match.group(5):
        name += '_concat'

    df_net = pd.read_csv(file, sep='\s+', lineterminator='\n')
    df_net = df_net.rename(columns={'pos': 'position', 'AA': 'residue', 'C': 'cleaved'})
    df_net['cleaved'] = df_net['cleaved'].replace(['S'], True).replace(['.'], False)

    peptides_from_netchop = calculate_cleaved_peptides.get_cleaved_peptides(df_net)

    intersec = [p for p in peptides_from_netchop if p in sars_peptides_netchop]
    intersection_count[name] = {'NetChop3.1': len(intersec)/len(peptides_from_netchop)}

for file in glob.glob('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/23AugEmbeddingPredictions/MHCI/pepsickle_17Aug_mhcI_*.txt'):
    match = re.search(r'_17Aug_(mhcI+(_\d+embedded)*(_hogvax)*)(_with_mutations)*', file)
    name = match.group(1).replace('_hogvax', '')

    df_pep = pd.read_csv(file, sep='\t', lineterminator='\n')
    df_pep = df_pep[df_pep['protein_id'] != 'HOGVAX']
    peptides_from_pepsickle = calculate_cleaved_peptides.get_cleaved_peptides(df_pep)

    intersec = [p for p in peptides_from_pepsickle if p in sars_peptides_pepsickle]
    if not name in intersection_count:
        print('ATTENTION! Is not in dict:', name)
        exit(-1)

    intersection_count[name]['Pepsickle'] = len(intersec) / len(peptides_from_pepsickle)

# add concat data
with open('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/23AugEmbeddingPredictions/MHCI/pepsickle_17Aug_mhcI_with_mutations.txt') as f:
    df_pep = pd.read_csv(file, sep='\t', lineterminator='\n')
    df_pep = df_pep[df_pep['protein_id'] != 'Concat']
    peptides_from_pepsickle = calculate_cleaved_peptides.get_cleaved_peptides(df_pep)

    intersec = [p for p in peptides_from_pepsickle if p in sars_peptides_pepsickle]
    if not name in intersection_count:
        print('ATTENTION! Is not in dict:', name)
        exit(-1)

    intersection_count['mhcI_concat']['Pepsickle'] = len(intersec) / len(peptides_from_pepsickle)


print(intersection_count)
df = pd.DataFrame(transform_data(intersection_count))
print(df)

# sorting = ['mhcII', 'mhcII_concat', 'mhcII_1embedded', 'mhcII_2embedded', 'mhcII_3embedded', 'mhcII_4embedded', 'mhcII_5embedded', 'mhcII_6embedded', 'mhcII_7embedded', 'mhcII_8embedded', 'mhcII_9embedded', 'mhcII_10embedded']
sorting = ['mhcI', 'mhcI_concat', 'mhcI_1embedded', 'mhcI_2embedded', 'mhcI_3embedded', 'mhcI_4embedded', 'mhcI_5embedded', 'mhcI_6embedded', 'mhcI_7embedded', 'mhcI_8embedded', 'mhcI_9embedded', 'mhcI_10embedded']
colors = ['#FFAD00', '#3A7A9B']
bar = alt.Chart(df).mark_bar(size=25).encode(
            x=alt.X('tool:O', axis=alt.Axis(title=None, labels=False, ticks=False)),
            y=alt.Y('fraction:Q', title='Fraction of cleaved peptides'),
            color=alt.Color('tool:O', scale=alt.Scale(range=colors), legend=alt.Legend(title='Tool')),
            column=alt.Column('embedding:O', title='Cleaved peptides that are equally predicted from SARS-CoV-2 proteome for different embedding lengths', sort=sorting)
        ).configure_axis(
            # grid=False
        ).configure_view(
            strokeOpacity=0
        ).properties(
            width=90
        )

bar.save('cleavage_plots/' + out_name + '.html')


# all_epi_ft = pd.read_pickle('/home/sara/Documents/VaccinesProject/ivp/Gifford_Data/AllEpitopeFeatures.pkl')
# known_peptides_pepsickle = [p for p in sars_peptides_pepsickle if p in all_epi_ft.index]
# print(len(known_peptides_pepsickle))
# known_peptides_netchop = [p for p in sars_peptides_netchop if p in all_epi_ft.index]
# print(len(known_peptides_netchop))

