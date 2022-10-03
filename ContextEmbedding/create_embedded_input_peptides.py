import pandas as pd
from ProcessData import read_peptides

max_embed = 10
mhc = 'mhcII'
peptide_range = [13, 26] # exclusive end value
original_feature_file = '../../Gifford_Data/AllEpitopeFeatures.pkl'
ba_file = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Data/25June_mhc2_netmhcii-4.1_pred_affinity_pivot_v1v2.pkl.gz'

for i in range(1, max_embed+1):
    embed = str(i)
    feature_file = '../../Gifford_Data/peptide_context_embedding/' + embed + 'EmbeddedEpitopeFeatures.pkl'
    out_file_ori = '../../Gifford_Data/peptides/' + mhc + '_peptides_with_mutations.pep'
    out_file = '../../Gifford_Data/peptides/' + mhc + '_peptides_with_mutations_' + embed + 'embedding.pep'

    ori_features = pd.read_pickle(original_feature_file)
    features = pd.read_pickle(feature_file)

    # filtering like Gifford
    selfp = open('../../Gifford_Data/self_pept.txt').read().split('\n')
    ori_peptides = ori_features[(ori_features['epi_len'].isin(range(*peptide_range))) &
                                (ori_features['glyco_probs']<=0.0) &
                                (ori_features['crosses_cleavage']==0) &
                                (ori_features['Is_self_pept']==False)].index.to_list()
    print(len(ori_peptides))

    # drop peptides that do not occur in ba predictions
    ba = pd.read_pickle(ba_file)
    na = []
    for pep in ori_peptides:
        if pep not in ba.index:
            ori_peptides.remove(pep)
            na.append(pep)
    print(len(na), na)
    print(len(ori_peptides))

    peptides = features[features['peptide'].isin(ori_peptides)].index
    print(len(peptides))

    with open(out_file, 'w') as ori_out:
        ori_out.write(('\n').join(peptides))

with open(out_file_ori, 'w') as out:
    out.write('\n'.join(ori_peptides))