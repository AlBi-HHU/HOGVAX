import pandas as pd
import re


def get_embedded_epitope_features(embed, all_epi_feat, min_len):
    with open('../../Gifford_Data/peptide_context_embedding/base_proteins_with_end_codon.fasta', 'r') as fasta:
        proteins = {}
        key = ''
        for line in fasta:
            if line.startswith('>'):
                key = re.match('> Protein: ([A-Z]+\d*[a-z]*) |', line).group(1)
                continue
            proteins[key] = line.strip('\n').strip('*')

    df_list = []
    for peptide in all_epi_feat.index:
        protein = all_epi_feat.loc[peptide]['protein']
        if protein == 'S1' or protein == 'S2':
            protein = 'S'
        pos = proteins[protein].index(peptide)
        start_embed = pos - embed
        end_embed = pos + len(peptide) + embed
        if start_embed < 0 or end_embed > len(proteins[protein]):
            if len(peptide) < min_len:
                embedded = peptide
            else:
                continue
        else:
            embedded = proteins[protein][start_embed:end_embed]
        # embedded = proteins[protein][max(pos-embed, 0):min(pos+len(peptide)+embed, len(proteins[protein]))]
        tmp = pd.concat([pd.Series({'embedded': embedded, 'peptide': peptide, 'embedding_length': embed}), all_epi_feat.loc[peptide]])
        df_list.append(tmp)

    new_df = pd.concat(df_list, axis=1).T.set_index('embedded')
    print('Writing ../../Gifford_Data/' + str(embed) + 'EmbeddedEpitopeFeatures.pkl')
    new_df.to_pickle('../../Gifford_Data/peptide_context_embedding/' + str(embed) + 'EmbeddedEpitopeFeatures.pkl')


max_embed = 10
min_len = 10
# read in all epitope features of original peptides
all_epi_feat = pd.read_pickle('../../Gifford_Data/AllEpitopeFeatures.pkl')

for i in range(1, max_embed+1):
    get_embedded_epitope_features(embed=i, all_epi_feat=all_epi_feat, min_len=min_len)