import gzip
import numpy as np
import pandas as pd
import find_superstrings
import read_peptides


def compare_row_entries(df, row_super):
    data = np.vstack([df.values, row_super.values])
    max_row = np.max(data, axis=0).reshape(row_super.shape)
    new_row = pd.DataFrame(max_row, index=row_super.index, columns=row_super.columns)
    return new_row


ba_file = snakemake.input['ba_predictions']
ba_df = pd.read_pickle(gzip.open(ba_file))
peptides = snakemake.params['pep']
print(len(peptides))

# first read in all peptides and identify superstrings
superstrings = find_superstrings.main(peptides)
print(len(superstrings))

ba_df.replace([np.inf, -np.inf], np.nan, inplace=True)
if ba_df.isnull().values.any():
    ba_df.fillna(0, inplace=True)
for super in superstrings:
    new_row = compare_row_entries(ba_df.loc[superstrings[super]], ba_df.loc[[super]])
    ba_df.loc[[super]] = new_row

ba_df.to_pickle(snakemake.output[0], compression='gzip')
