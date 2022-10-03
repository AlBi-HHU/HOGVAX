import os
import pandas as pd
import numpy as np


def embedding_ba_predictions(embedded_epitope_features, ba_affinities, embedded_peptides, path):
    all_embedded = pd.read_pickle(embedded_epitope_features)

    mapped_peptides = pd.DataFrame(all_embedded.loc[embedded_peptides]['peptide']).reset_index(level=0).set_index('peptide')
    same_pep = list(set([p for p in mapped_peptides.index if p in ba_affinities.index]))
    mapped_peptides = mapped_peptides.loc[same_pep]
    filtered = ba_affinities.reset_index()
    filtered.Peptide = filtered.Peptide.apply(lambda x: mapped_peptides.loc[x]['embedded'] if x in same_pep else np.nan)
    filtered = filtered[filtered['Peptide'].notna()]
    filtered = filtered.set_index(['Peptide'])
    for embed in mapped_peptides['embedded']:
        if embed not in filtered.index:
            print(embed, ' is not in filtered index!!!')
            exit(-1)

    if not os.path.exists(path + 'data/'):
        os.mkdir(path + 'data/')
    filtered.to_pickle(path + 'data/' + str(len(embedded_peptides)) + '_embedding_ba_predictions.pkl.gz', compression='gzip')
    return filtered
