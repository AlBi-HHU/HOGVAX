import sys
import pandas as pd


def get_cleaved_peptides(cleavage):
    vax_seq = ''.join(cleavage['residue'])
    start = 0
    peptides_from_prediction = []
    # print(cleavage.iloc[0])
    # print(cleavage[cleavage['cleaved']])
    for outer in cleavage[cleavage['cleaved']].itertuples():
        outer = outer._asdict()
        # cut happens after cleaved = True but file count is 1-based
        next_start = outer['position']
        for inner in cleavage[(cleavage['cleaved']) & (cleavage['position'] > start)].itertuples():
            inner = inner._asdict()
            # 1-based!
            pos = inner['position']
            peptide = vax_seq[start:pos]
            if 8 <= len(peptide) <= 25:
                peptides_from_prediction.append(vax_seq[start:pos])
        start = next_start

    return list(set(peptides_from_prediction))


def main(pepsickle_file, netchop_file):
    with open(pepsickle_file, 'r') as file:
        df_pep = pd.read_csv(file, sep='\t', lineterminator='\n')

    peptides_from_pepsickle = get_cleaved_peptides(df_pep)

    with open('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/SarsCoVCleavage/pepsickle_sars_cov_2_cleavage_peptides.txt', 'w') as out:
        out.write('\n'.join(peptides_from_pepsickle))

    with open(netchop_file, 'r') as file:
        df_net = pd.read_csv(file, sep='\s+', lineterminator='\n')

    df_net = df_net.rename(columns={'pos': 'position', 'AA': 'residue', 'C': 'cleaved'})
    df_net['cleaved'] = df_net['cleaved'].replace(['S'], True).replace(['.'], False)

    peptides_from_netchop = get_cleaved_peptides(df_net)

    with open(
            '/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/SarsCoVCleavage/netchop_sars_cov_2_cleavage_peptides.txt',
            'w') as out:
        out.write('\n'.join(peptides_from_netchop))


if __name__ == '__main__':
    pepsickle = sys.argv[1]
    netchop = sys.argv[2]
    main(pepsickle, netchop)

