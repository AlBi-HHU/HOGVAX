import pandas as pd
import numpy as np
from ProcessData import read_peptides

peptide_file = '../../Gifford_Data/peptides/mhcI_peptides_with_mutations_1embedding.pep'
embedded_epitope_features = '../../Gifford_Data/peptide_context_embedding/1EmbeddedEpitopeFeatures.pkl'
ba_file = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Data/25June_mhc1_netmhc-4.1_pred_affinity_pivot.pkl.gz'

embedded_peptides = read_peptides.main(peptide_file)
# embedded_peptides = ['MYSFVSEET', 'MYSFVSEETG', 'YSFVSEETGT', 'SFVSEETGTL']
all_embedded = pd.read_pickle(embedded_epitope_features)
ba_affinities = pd.read_pickle(ba_file)

mapped_peptides = pd.DataFrame(all_embedded.loc[embedded_peptides]['peptide']).reset_index(level=0).set_index('peptide')
print(mapped_peptides)
same_pep = list(set([p for p in mapped_peptides.index if p in ba_affinities.index]))
print(same_pep)
mapped_peptides = mapped_peptides.loc[same_pep]
filtered = ba_affinities.reset_index()
filtered.Peptide = filtered.Peptide.apply(lambda x: mapped_peptides.loc[x]['embedded'] if any(x == p for p in same_pep) else np.nan)
filtered = filtered[filtered['Peptide'].notna()]
print(filtered)
filtered = filtered.set_index(['Peptide'])
for embed in mapped_peptides['embedded']:
    if embed not in filtered.index:
        print(embed, ' is not in filtered index!!!')
        exit(-1)

print(filtered)
# new_df.to_pickle('../../Gifford_Data/' + str(embed) + 'embedded_25June_mhc1_netmhc-4.1_pred_affinity_pivot.pkl.gz', compression='gzip')