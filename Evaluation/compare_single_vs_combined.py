# used this script for debugging
import gzip
import pandas as pd
from ProcessData import read_peptides, binarize_entries, find_superstrings

combined = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/13July_mhcIandII/ilp/pep_out/41947_chosen_peptides_hog.txt'
combined_substr = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/13July_mhcIandII/ilp/pep_out/41947_chosen_peptides_hog_inc_substrings.txt'
combined_peptides = read_peptides.main(combined)
combined_peptides_substr = read_peptides.main(combined_substr)
print(len(combined_peptides))

# mhcI = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/11July_filt_substrmhcI/ilp/pep_out/4512_chosen_peptides_hog_inc_substrings.txt'
# mhcI_peptides = read_peptides.main(mhcI)
# mhcII = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/12July_filt_substrmhcII/ilp/pep_out/37435_chosen_peptides_hog_inc_substrings.txt'
# mhcII_peptides = read_peptides.main(mhcII)

# mhcI_equal = set(mhcI_peptides).intersection(combined_peptides)
# print(len(mhcI_equal), 'of', len(mhcI_peptides))
#
# mhcII_equal = set(mhcII_peptides).intersection(combined_peptides)
# print(len(mhcII_equal), 'of', len(mhcII_peptides))

combined_ba = pd.read_pickle(gzip.open('/home/sara/Documents/VaccinesProject/ivp/Gifford_Data/25June_combined_mhc_netmhc_pred_affinity.pkl.gz'))
combined_ba_substr = pd.read_pickle(gzip.open('/home/sara/Documents/VaccinesProject/ivp/Gifford_Data/substr_modified_25June_combined_mhc_netmhc_pred_affinity.pkl.gz'))
# mhcI_ba = pd.read_pickle(gzip.open('../../Gifford_Data/25June_mhc1_netmhc-4.1_pred_affinity_pivot.pkl.gz'))
# mhcII_ba = pd.read_pickle(gzip.open('../../Gifford_Data/25June_mhc2_netmhcii-4.1_pred_affinity_pivot_v1v2.pkl.gz'))

combined_af = pd.read_pickle('../../Gifford_Data/IEDB_combined_population_frequency.pkl')
# mhcI_af = pd.read_pickle('../../Gifford_Data/IEDB_population_frequency2392_normalized.pkl')
# mhcII_af = pd.read_pickle('../../Gifford_Data/IEDB_population_frequency_mhc2_275normalized.pkl')

print(set(combined_peptides_substr).intersection(combined_peptides))

bin_ba = binarize_entries.binarize_entries(combined_ba, 0.638)
bin_ba = bin_ba.loc[combined_peptides_substr]
hit_alleles = bin_ba[[x for x in bin_ba if sum(bin_ba[x]) > 0]].columns
score = sum(combined_af[[al for al in hit_alleles if al in combined_af]].loc['World'])
print(score)

bin_ba_sub = binarize_entries.binarize_entries(combined_ba_substr, 0.638)
bin_ba_sub = bin_ba_sub.loc[combined_peptides]
hit_alleles_sub = bin_ba_sub[[x for x in bin_ba_sub if sum(bin_ba_sub[x]) > 0]].columns
score_sub = sum(combined_af[[al for al in hit_alleles_sub if al in combined_af]].loc['World'])
print(score_sub)

all_peptides = read_peptides.main('/home/sara/Documents/VaccinesProject/ivp/Gifford_Data/peptides/combined_filtered_mhc_peptides_gifford_style.pep')
# superstrings = find_superstrings.main(all_peptides)
super_covered_alleles = []
for al in hit_alleles_sub:
    if al not in hit_alleles:
        super_covered_alleles.append(al)
        # for p in bin_ba_sub.index:
        #     if bin_ba_sub[al][p] > 0.0 and combined_ba[al][p] == 0.0:
        #         print('Peptide {} covering allele {} in merged data, but not in original data'.format(p, al))
        #         substrings = superstrings[p]
        #         print(substrings)
        #         if not any(combined_ba[al][x] > 0.0 for x in substrings):
        #             print('None of the substrings has non-zero value in original data')
        #         else:
        #             for x in substrings:
        #                 if x not in combined_peptides_substr:
        #                     print('why is {} not in combined substr peptied, but {} is in and superstr'.format(x, p))

print(super_covered_alleles)