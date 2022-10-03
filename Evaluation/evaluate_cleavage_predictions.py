import sys
import pandas as pd
from ProcessData import read_peptides
from EvaluateVaccines import calculate_cleaved_peptides


def count_pred_peptides(chosen_peptides, peptides_from_prediction):
    counter = 0
    peptides = []
    for pep in chosen_peptides:
        if pep in peptides_from_prediction:
            counter += 1
            peptides.append(pep)

    print(len(chosen_peptides), counter, sorted(peptides))
    return counter, peptides


# read peptides chosen by ILP including all substrings
peptide_file = sys.argv[1]
chosen_peptides = read_peptides.main(peptide_file)

# Note that netchop requires manually deleting header and footer lines + last column modification
pepsickle_file = sys.argv[2]
netchop_file = sys.argv[3]

# csv file to write output
out_file = sys.argv[4]

concat_eval = False
if len(sys.argv) > 5:
    concat_eval = True

with open(pepsickle_file, 'r') as file:
    df_pep = pd.read_csv(file, sep='\t', lineterminator='\n')
    if concat_eval:
        df_pep = df_pep[df_pep['protein_id'] != 'HOGVAX']
    else:
        df_pep = df_pep[df_pep['protein_id'] != 'Concat']

peptides_from_pepsickle = calculate_cleaved_peptides.get_cleaved_peptides(df_pep)
count_pepsickle, pred_peptides_pepsickle = count_pred_peptides(chosen_peptides, peptides_from_pepsickle)


with open(netchop_file, 'r') as file:
    df_net = pd.read_csv(file, sep='\s+', lineterminator='\n')

df_net = df_net.rename(columns={'pos': 'position', 'AA': 'residue', 'C': 'cleaved'})
df_net['cleaved'] = df_net['cleaved'].replace(['S'], True).replace(['.'], False)

peptides_from_netchop = calculate_cleaved_peptides.get_cleaved_peptides(df_net)
count_netchop, pred_peptides_netchop = count_pred_peptides(chosen_peptides, peptides_from_netchop)

d = {'Chosen Peptides': chosen_peptides,
     'Pepsickle Predicted Peptides': pred_peptides_pepsickle,
     'NetChop Predicted Peptides': pred_peptides_netchop}
print(d)
df = pd.DataFrame.from_dict(d, orient='index')
df = df.T
df.to_csv(out_file)



