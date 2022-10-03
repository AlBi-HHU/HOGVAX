import glob
import re
import pandas as pd
import altair as alt

out_name = '24Sept_mhcI_cleavage_prediction_for_embedding_methods'

data = {'embedding': [], 'tool': [], 'fraction': []}
for file in glob.glob('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/24SeptEmbeddingPredictions/MHCI/CleavagePredictionsForEvaluation/*.csv'):
    match = re.search(r'_17Aug_(mhcI+(_\d+embedded)*_+hogvax)(_with_mutations)*(_concat)*', file)
    name = match.group(1).replace('_hogvax', '')
    if match.group(4):
        name += '_concat'

    tmp = pd.read_csv(file)
    tmp = tmp.drop(columns=tmp.columns[0])
    chosen_pep = len(tmp['Chosen Peptides'].dropna())
    pepsickle_pep = len(tmp['Pepsickle Predicted Peptides'].dropna())
    netchop_pep = len(tmp['NetChop Predicted Peptides'].dropna())

    check_a = [p for p in tmp['Pepsickle Predicted Peptides'] if p in list(tmp['Chosen Peptides'])]
    if len(check_a) != pepsickle_pep:
        print('Attention!')
        exit(-1)

    check_b = [p for p in tmp['NetChop Predicted Peptides'] if p in list(tmp['Chosen Peptides'])]
    if len(check_b) != netchop_pep:
        print('Attention!')
        exit(-1)

    data['embedding'].extend([name, name])
    data['tool'].extend(['NetChop3.1', 'Pepsickle'])
    data['fraction'].extend([netchop_pep/chosen_pep, pepsickle_pep/chosen_pep])

data = pd.DataFrame.from_dict(data)
print(data)

# sorting = ['mhcII', 'mhcII_concat', 'mhcII_1embedded', 'mhcII_2embedded', 'mhcII_3embedded', 'mhcII_4embedded', 'mhcII_5embedded', 'mhcII_6embedded', 'mhcII_7embedded', 'mhcII_8embedded', 'mhcII_9embedded', 'mhcII_10embedded']
sorting = ['mhcI', 'mhcI_concat', 'mhcI_1embedded', 'mhcI_2embedded', 'mhcI_3embedded', 'mhcI_4embedded', 'mhcI_5embedded', 'mhcI_6embedded', 'mhcI_7embedded', 'mhcI_8embedded', 'mhcI_9embedded', 'mhcI_10embedded']
colors = ['#DA0034', '#3A7A9B']
bar = alt.Chart(data).mark_bar(size=25).encode(
            x=alt.X('tool:O', axis=alt.Axis(title=None, labels=False, ticks=False)),
            y=alt.Y('fraction:Q', title='Cleaved peptides'),
            color=alt.Color('tool:O', scale=alt.Scale(range=colors), legend=alt.Legend(title='Tool')),
            column=alt.Column('embedding:O', title='Cleaved peptides for different context embedding lengths', sort=sorting)
        ).configure_axis(
            # grid=False
        ).configure_view(
            strokeOpacity=0
        ).properties(
            width=90
        )

bar.save('cleavage_plots/' + out_name + '.html')