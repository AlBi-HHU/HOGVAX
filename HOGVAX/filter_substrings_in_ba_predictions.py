import os
import numpy as np
import pandas as pd


def modify_ba_predictions(superstrings, ba_df, pep_count, path):
    ba_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    if ba_df.isnull().values.any():
        ba_df.fillna(0, inplace=True)
    # only works for binary values
    for super in superstrings:
        for sub in superstrings[super]:
            ba_df.loc[super] = ba_df.loc[super] | ba_df.loc[sub]

    if not os.path.exists(path + 'data/'):
        os.mkdir(path + 'data/')
    ba_df.to_pickle(path + 'data/' + str(pep_count) + '_substr_modified_ba_predictions.pkl.gz', compression='gzip')

    return ba_df
