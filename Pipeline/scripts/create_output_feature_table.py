import sys
import pandas as pd


def collect_features(peptides, all_features):
    ilp_out_features = []
    for pep in peptides:
        all = all_features.loc[pep]
        dict = {'Peptide': pep, 'Protein': all['protein'],
                'Start_pos': all['start_pos'], 'Pep_len': len(pep),
                'Mutation_prob': all['perc_mutated'], 'Glyco_prob': all['glyco_probs'],
                'Cleavage_prob': all['crosses_cleavage'], 'SARS-CoV-1_peptide': all['In_SARS_Cov1']}
        ilp_out_features.append(dict)
    return ilp_out_features


def read_peptides(file):
    peptides = []
    with open(file) as f:
        for line in f:
            peptides.append(line.strip('\n'))
    return peptides


# all epitope freatures based on gifford data
all_epi_features_file = sys.argv[1]
# ilp peptide output
peptides_mhcI_file = sys.argv[2]
peptides_mhcII_file = sys.argv[3]

all_epi_features = pd.read_pickle(all_epi_features_file)
peptides_mhcI = read_peptides(peptides_mhcI_file)
peptides_mhcII = read_peptides(peptides_mhcII_file)

with pd.ExcelWriter('peptides_output/all_peptide_vaccine_features.xlsx', engine='openpyxl') as writer:
    features = collect_features(peptides_mhcI, all_epi_features)
    df = pd.DataFrame(features)
    df.to_excel(writer, sheet_name='Sheet1')

    features = collect_features(peptides_mhcII, all_epi_features)
    df = pd.DataFrame(features)
    df.to_excel(writer, sheet_name='Sheet2')
