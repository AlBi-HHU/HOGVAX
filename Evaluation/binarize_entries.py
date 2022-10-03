import pandas as pd


def binarize_entries(df, th):
    for key in df.keys():
        df.loc[df[key] >= th, key] = 1.0
        df.loc[df[key] < th, key] = 0.0
    return df


def set_small_to_zero(df, th):
    for key in df.keys():
        df.loc[df[key] < th, key] = 0.0
    return df
