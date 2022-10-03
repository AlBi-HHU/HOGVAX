import os
import numpy as np
import pandas as pd
import find_superstrings


def compare_row_entries(df, row_super):
    data = np.vstack([df.values, row_super.values])
    max_row = np.max(data, axis=0).reshape(row_super.shape)
    new_row = pd.DataFrame(max_row, index=row_super.index, columns=row_super.columns)
    return new_row


def modify_ba_predictions(peptides, ba_df, path):
    # first read in all peptides and identify superstrings
    superstrings = find_superstrings.find_superstrings(peptides)
    print(len(superstrings))

    ba_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    if ba_df.isnull().values.any():
        ba_df.fillna(0, inplace=True)
    for super in superstrings:
        new_row = compare_row_entries(ba_df.loc[superstrings[super]], ba_df.loc[[super]])
        ba_df.loc[[super]] = new_row

    # print(superstrings)
    if not os.path.exists(path + 'data/'):
        os.mkdir(path + 'data/')
    ba_df.to_pickle(path + 'data/' + str(len(peptides)) + '_substr_modified_ba_predictions.pkl.gz', compression='gzip')

    return ba_df
