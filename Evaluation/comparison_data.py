import gzip
import pickle
import os
import pandas as pd
import numpy as np
import evaluate_vaccine as ev
import plot_evaluation_results as plot
import read_peptides, binarize_entries, create_haplotype_ba_predictions
from collections import Counter


def get_cmp_obj(peptides, all_ba, af, smaller_to_zero=False, binary=False, filter=True):
    # only keep peptide columns corresponding to chosen peptides that also exist in prediction file
    good_keys = all_ba.index.intersection(peptides)
    pep_ba = all_ba.loc[good_keys]

    if smaller_to_zero:
        pep_ba = binarize_entries.set_small_to_zero(pep_ba, 0.638)

    if binary:
        pep_ba = binarize_entries.binarize_entries(pep_ba, 0.638)

    if filter:
        if len(pep_ba.columns) != len(af.columns):
            columns_to_keep = [key for key in pep_ba if key in af]
            pep_ba = pep_ba[columns_to_keep]
            af = af[columns_to_keep]

    return ev.EvaluationMetrices(peptides, pep_ba, af)
    

def collect_data_for_bar_plot(method_A, method_A_mhc1, method_A_mhc2, method_B, method_B_mhc1, method_B_mhc2,
                              mhc1_all_ba, mhc1_af, mhc2_all_ba, mhc2_af, outdir):

    results = {'Method': [], 'Population Coverage': [], 'Number of Peptides': []}

    # mhc1
    peptides = read_peptides.main(method_A_mhc1)
    results['Method'].append('MHCI ' + method_A)
    results['Population Coverage'].append(get_cmp_obj(peptides, mhc1_all_ba, mhc1_af, smaller_to_zero=True).pred_at_least_one('World'))
    results['Number of Peptides'].append(len(peptides))

    peptides = read_peptides.main(method_B_mhc1)
    results['Method'].append('MHCI ' + method_B)
    results['Population Coverage'].append(get_cmp_obj(peptides, mhc1_all_ba, mhc1_af, smaller_to_zero=True).pred_at_least_one('World'))
    results['Number of Peptides'].append(len(peptides))

    # mhc2
    peptides = read_peptides.main(method_A_mhc2)
    results['Method'].append('MHCII ' + method_A)
    results['Population Coverage'].append(get_cmp_obj(peptides, mhc2_all_ba, mhc2_af, smaller_to_zero=True).pred_at_least_one('World'))
    results['Number of Peptides'].append(len(peptides))

    peptides = read_peptides.main(method_B_mhc2)
    results['Method'].append('MHCII ' + method_B)
    results['Population Coverage'].append(get_cmp_obj(peptides, mhc2_all_ba, mhc2_af, smaller_to_zero=True).pred_at_least_one('World'))
    results['Number of Peptides'].append(len(peptides))

    df = pd.DataFrame.from_dict(results)
    plot.EvaluationPlots(df, outdir).ilp_opitvax_bar_plot()
    df.to_csv(outdir + 'evalvax-unlinked.csv')


def collect_data_for_locus_cov_pie(peptide_file, mhc, mhc_ba, af, loci, population, method, outdir):
    cov_loci = {}
    hit_loci = {}
    hit_loci_freq = {}
    peptides = read_peptides.main(peptide_file)
    for locus in loci:
        covered, unpred_alleles, unpred_frac, total, allele_hit_counts, hit_allele_freq = get_cmp_obj(peptides,
                                                                                                      mhc_ba,
                                                                                                      af,
                                                                                                      binary=True,
                                                                                                      filter=False
                                                                                                      ).number_covered_alleles(population, locus)
        uncovered = abs(total - covered - unpred_frac)
        values = [round(v, 4) for v in [covered, unpred_frac, uncovered]]
        cov_loci[locus] = {'allele_category': ['Covered', 'No-BA-prediction', 'Uncovered'], 'values': values}
        hit_loci[locus] = allele_hit_counts

        categories = []
        fractions = []
        for key in sorted(hit_allele_freq.keys(), key=float):
            categories.append(str(key))
            fractions.append(hit_allele_freq[key])
        hit_loci_freq[locus] = {'allele_category': categories + ['No-BA-prediction', 'Uncovered'],
                                'values': fractions + [unpred_frac, uncovered],
                                'method': [locus + method] * (len(categories) + 2)}
        print(hit_loci_freq)

    plot.EvaluationPlots(cov_loci, outdir).covered_locus_pie(loci, mhc, method)
    plot.EvaluationPlots(hit_loci, outdir).hit_pies(loci, mhc, method)
    plot.EvaluationPlots(hit_loci_freq, outdir).hit_covered_locus_pie(loci, mhc, method)
    return hit_loci_freq


def collect_data_for_haplotype_eval(peptide_file, single_ba, hap_ba, af, max_int, mhc, method, outdir):
    peptides = read_peptides.main(peptide_file)
    single_ba = binarize_entries.binarize_entries(single_ba, 0.638).loc[[p for p in peptides if p in single_ba.index]]
    hit_dict = {}
    hit_n_dict = []
    for population in af.index:
        # hits is number of hits per genotype, hits_n is probability for n>=N hits per person in population
        hits, hits_n, exp, no_predictions = get_cmp_obj(peptides, hap_ba.T, af, binary=True).pred_genotype_hits(population, single_ba, max_int, mhc)
        hit_dict[population] = [hits, {'expected': [exp]}]
        tmp = pd.DataFrame(hits)
        tmp.to_csv(outdir + 'hit_hist_' + mhc + '_' + method + '.csv')
        # preparation for grouped bar plot
        hit_n_dict.append(hits_n | {'population': [population]*max_int})
    # plot histogram of peptide-hla hits
    plot.EvaluationPlots(hit_dict, outdir).hit_histogram(mhc, method)

    cur_probs = [np.array(dic['prob']) for dic in hit_n_dict]
    mean_probs = list(np.mean(cur_probs, axis=0))
    hit_n_dict.append({'prob': mean_probs, 'n': [i for i in range(1, max_int+1)], 'population': ['Average']*max_int})

    # reformat dict for grouped bar plots
    df = pd.concat(pd.DataFrame(entry) for entry in hit_n_dict)
    plot.EvaluationPlots(df, outdir).population_hit_bar_plot(mhc, method)
    df.to_csv(outdir + 'min_hit_dict' + mhc + '_' + method + '.csv')


def compute_hap_freq_no_ld(population, allele_freq, loci):
    hap_data_values = []
    hap_data_labels = []
    for locus in loci:
        tmp = allele_freq[locus].loc[population]
        tmp = tmp[tmp > 0]
        if 'unknown' in tmp:
            tmp = tmp.drop('unknown')
        hap_data_labels.append(tmp.index)
        hap_data_values.append(tmp.values)

    hap_freq = []
    for i in range(len(hap_data_labels[0])):
        print(i)
        for j in range(len(hap_data_labels[1])):
            for k in range(len(hap_data_labels[2])):
                f = hap_data_values[0][i] * hap_data_values[1][j] * hap_data_values[2][k]
                hap_freq.append(((hap_data_labels[0][i], hap_data_labels[1][j], hap_data_labels[2][k]), f))

    unknown = 1 - sum(i for _, i in hap_freq)
    hap_freq.append((('unknown', 'unknown', 'unknown'), unknown))
    df = pd.DataFrame(hap_freq, columns=['Haplotype', population])
    df = df.set_index('Haplotype').T

    return df


def collect_data_for_geographic_eval(threads, peptide_file, populations, allele_freq, single_ba, loci, mhc_type):
    peptides = read_peptides.main(peptide_file)
    hit_n_dict = []
    for population in populations:
        print(population)
        haplotype_data = compute_hap_freq_no_ld(population, allele_freq, loci)
        print('Calculated haplotype frequencies for ' + population)
        hap_ba = create_haplotype_ba_predictions.main(threads, single_ba, haplotype_data, loci, population, mhc_type)
        print('Calculated haplotype binding affinities')
        single_ba = binarize_entries.binarize_entries(single_ba, 0.638).loc[peptides]
        hits, hits_n, exp, no_predictions = get_cmp_obj(peptides, hap_ba.T, haplotype_data, binary=True).pred_genotype_hits(population, single_ba, 1, mhc_type)
        print('Calculated evaluation results')
        hit_n_dict.append(hits_n | {'population': population})

    with open('data/hit_predictions_no_ld.pkl', 'wb') as file:
        pickle.dump(hit_n_dict, file)


def collect_data_for_composition_eval(peptide_file, feature_table, method, mhc_type, outdir):
    peptides = read_peptides.main(peptide_file)
    # composition = dict(Counter(feature_table.loc[peptides].protein))
    composition = dict(Counter(feature_table.loc[peptides].protein))
    df = pd.DataFrame(composition.items(), columns=['Protein', 'Count'])
    print(df)
    plot.EvaluationPlots(df, outdir).composition_plot(mhc_type, method)


def call_single_evaluations(peptides, mhc_type, feature_table, allele_ba_df, haplotype_ba_df, allele_af_df, haplotype_af_df, loci,
                            population, name, outdir, max_hit=10):
    # print('Create pie charts')
    # hit_dict = collect_data_for_locus_cov_pie(peptides, mhc_type, allele_ba_df, allele_af_df, loci, population, name, outdir)
    print('Create histograms for haplotype comparison')
    collect_data_for_haplotype_eval(peptides, allele_ba_df, haplotype_ba_df, haplotype_af_df, max_hit, mhc_type, name, outdir)
    # print('Plot vaccine composition')
    # collect_data_for_composition_eval(peptides, feature_table, mhc_type, name, outdir)
    # return hit_dict


def new_df(dict_list):
    new_dict = {'allele_category': [], 'values': [], 'method': []}
    for dict in dict_list:
        for key in dict:
            new_dict['allele_category'].extend(dict[key]['allele_category'])
            new_dict['values'].extend(dict[key]['values'])
            new_dict['method'].extend(dict[key]['method'])
    df = pd.DataFrame(new_dict)
    return df


def main():
    # optivax_mhc1 = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc1_optivax_unlinked_World_June28/28June_optivax_unlinked_mhc1_World_beam_peptides.txt'
    # optivax_mhc2 = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc2_optivax_unlinked_World_June28/28June_optivax_unlinked_mhc2_World_beam_peptides.txt'
    #
    # hogvax_mhc1 = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/18Sept_mhcI_122peptides/ilp/pep_out/4512_chosen_peptides_hog_inc_substrings.txt'
    # hogvax_mhc2 = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/18Sept_mhcII_1112/ilp/pep_out/37435_chosen_peptides_hog_inc_substrings.txt'

    # haplotypes
    optivax_mhc1 = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc1_optivax_robust_haplotype_modifications_July27/27July_optivax_robust_mhc1_beam_peptides.txt'
    optivax_mhc2 = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc2_optivax_robust_haplotype_modifications_July27/27July_optivax_robust_mhc2_beam_peptides.txt'

    hogvax_mhc1 = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/28July_mhcI_haplotypes/ilp/pep_out/4519_chosen_peptides_hog_inc_substrings.txt'
    hogvax_mhc2 = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/28July_mhcII_haplotypes/ilp/pep_out/37435_chosen_peptides_hog_inc_substrings.txt'

    # hogvax_mhc1_gen = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/29Aug_mhcI_100k_genotypes_average/ilp/pep_out/4519_chosen_peptides_hog_inc_substrings.txt'
    # hogvax_mhc2_gen = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/29Aug_mhcII_100k_genotypes_average/ilp/pep_out/37435_chosen_peptides_hog_inc_substrings.txt'
    #
    # optivax_mhc1_gen = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc1_optivax_robust_genotypes_July31/31July_optivax_robust_mhc1_genotype_peptides.txt'
    # optivax_mhc2_gen = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc2_optivax_robust_genotypes_July31/31July_optivax_robust_mhc2_genotype_peptides.txt'

    # optivax_mhc1_spike = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc1_optivax_unlinked_spike_protein/5Sep_optivax_unlinked_spike_peptides.txt'
    # optivax_mhc2_spike = '/home/sara/Documents/VaccinesProject/ivp/Gifford_Code/mhc2_optivax_unlinked_spike_protein/5Sep_optivax_unlinked_spike_peptides.txt'
    #
    # hogvax_mhc1_spike = '/home/sara/Documents/VaccinesProject/ivp/Code/HOGVAX/spike_new_mhc1/pep_out/639_chosen_peptides_hogvax_inc_substrings.txt'
    # hogvax_mhc2_spike = '/home/sara/Documents/VaccinesProject/ivp/Code/HOGVAX/spike_new_mhc2/pep_out/5064_chosen_peptides_hogvax_inc_substrings.txt'

    # concat_mhcI = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/28July_mhcI/new_concat(not using mod ba predictions)/naively_chosen_peptides.txt'
    # concat_mhcII = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/28July_mhcII/new_concat(no mod ba data)/naively_chosen_peptides.txt'

    # evaluation for combined mhcI and II vaccine
    # comb_peptide_file = '/home/sara/Documents/VaccinesProject/ivp/Pipeline/25Sept_mhcIandII/ilp/pep_out/41947_chosen_peptides_hog_inc_substrings.txt'

    all_epitope_features = pd.read_pickle('../../Gifford_Data/AllEpitopeFeatures.pkl')

    mhc1_all_ba = pd.read_pickle(gzip.open('../../Gifford_Data/25June_mhc1_netmhc-4.1_pred_affinity_pivot.pkl.gz'))
    mhc1_af = pd.read_pickle('../../Gifford_Data/IEDB_population_frequency2392_normalized.pkl')

    mhc2_all_ba = pd.read_pickle(gzip.open('../../Gifford_Data/25June_mhc2_netmhcii-4.1_pred_affinity_pivot_v1v2.pkl.gz'))
    mhc2_af = pd.read_pickle('../../Gifford_Data/IEDB_population_frequency_mhc2_275normalized.pkl')

    mhc1_haplotype_ba = pd.read_pickle(gzip.open('../../Gifford_Data/all_mhc1_predictions_for_haplotypes.pkl.gz')).T
    mhc1_haplotype_freq = pd.read_pickle('../../Gifford_Data/haplotype_frequency_marry.pkl')

    mhc2_haplotype_ba = pd.read_pickle(gzip.open('../../Gifford_Data/all_mhc2_predictions_for_haplotypes.pkl.gz')).T
    mhc2_haplotype_freq = pd.read_pickle('../../Gifford_Data/haplotype_frequency_marry2.pkl')

    outdir = 'plots_haplotypes/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # collect data and create bar plot
    print('Create simple bar plot')
    collect_data_for_bar_plot('OptiVax-Unlinked', optivax_mhc1, optivax_mhc2, 'HOGVAX', hogvax_mhc1, hogvax_mhc2,
                              mhc1_all_ba, mhc1_af, mhc2_all_ba, mhc2_af, outdir)
    # collect_data_for_bar_plot('Concat', concat_mhcI, concat_mhcII, 'HOGVAX', hogvax_mhc1, hogvax_mhc2,
    #                           mhc1_all_ba, mhc1_af, mhc2_all_ba, mhc2_af, outdir)
    # collect_data_for_bar_plot('Combined', comb_peptide_file, comb_peptide_file, 'HOGVAX', hogvax_mhc1, hogvax_mhc2,
    #                           mhc1_all_ba, mhc1_af, mhc2_all_ba, mhc2_af, outdir)

    print('Call eval for MHC1 HOGVAX')
    hits_alleles_hogvax_mhc1 = call_single_evaluations(hogvax_mhc1, 'mhc1', all_epitope_features,
                                                       mhc1_all_ba, mhc1_haplotype_ba, mhc1_af, mhc1_haplotype_freq,
                                                       ['HLA-A', 'HLA-B', 'HLA-C'], 'World', 'HOGVAX', outdir)
    print('Call eval for MHC1 OptiVax')
    hits_alleles_optivax_mhc1 = call_single_evaluations(optivax_mhc1, 'mhc1', all_epitope_features,
                                                        mhc1_all_ba, mhc1_haplotype_ba, mhc1_af, mhc1_haplotype_freq,
                                                        ['HLA-A', 'HLA-B', 'HLA-C'], 'World', 'OptiVax', outdir)
    print('Call eval for MHC2 HOGVAX')
    hits_alleles_hogvax_mhc2 = call_single_evaluations(hogvax_mhc2, 'mhc2', all_epitope_features,
                                                       mhc2_all_ba, mhc2_haplotype_ba, mhc2_af, mhc2_haplotype_freq,
                                                       ['DRB1', 'HLA-DP', 'HLA-DQ'], 'World', 'HOGVAX', outdir)
    print('Call eval for MHC2 OptiVax')
    hits_alleles_optivax_mhc2 = call_single_evaluations(optivax_mhc2, 'mhc2', all_epitope_features,
                                                        mhc2_all_ba, mhc2_haplotype_ba, mhc2_af, mhc2_haplotype_freq,
                                                        ['DRB1', 'HLA-DP', 'HLA-DQ'], 'World', 'OptiVax', outdir)

    # print('Plot Allele hit frequencies hogvax + optivax mhc2')
    # df = new_df([hits_alleles_hogvax_mhc2, hits_alleles_optivax_mhc2])
    # plot.EvaluationPlots(df, outdir).stacked_hit_covered_alleles('mhc2', 'hogvax_optivax_unlinked')

    # call_single_evaluations(concat_mhcI, 'mhc1', all_epitope_features, mhc1_all_ba, mhc1_haplotype_ba, mhc1_af, mhc1_haplotype_freq,
    #                         ['HLA-A', 'HLA-B', 'HLA-C'], 'World', 'Concat')
    # call_single_evaluations(concat_mhcII, 'mhc2', all_epitope_features, mhc2_all_ba, mhc2_haplotype_ba, mhc2_af, mhc2_haplotype_freq,
    #                         ['DRB1', 'HLA-DP', 'HLA-DQ'], 'World', 'Concat')
    #
    # call_single_evaluations(comb_peptide_file, 'mhc1', all_epitope_features, mhc1_all_ba, mhc1_haplotype_ba, mhc1_af,
    #                         mhc1_haplotype_freq, ['HLA-A', 'HLA-B', 'HLA-C'], 'World', 'Combined')
    # call_single_evaluations(comb_peptide_file, 'mhc2', all_epitope_features, mhc2_all_ba, mhc2_haplotype_ba, mhc2_af,
    #                         mhc2_haplotype_freq, ['DRB1', 'HLA-DP', 'HLA-DQ'], 'World', 'Combined')
    # countries = ['Europe',
    #              'Northeast Asia',
    #              'Oceania', 'East Africa',
    #              'East Asia', 'South Africa',
    #              'Southeast Asia',
    #              'Central Africa', 'North America',
    #              'South America', 'West Africa',
    #              'West Indies', 'South Asia',
    #              'North Africa', 'Southwest Asia']
    #
    # threads = 1
    # collect_data_for_geographic_eval(threads, ilp_mhc1, countries, mhc1_af, mhc1_all_ba, ['HLA-A', 'HLA-B', 'HLA-C'], 'mhc1')

if __name__ == '__main__':
    main()
