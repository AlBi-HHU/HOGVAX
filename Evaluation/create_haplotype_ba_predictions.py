import argparse
import sys
import pandas as pd
from ProcessData import read_peptides
from multiprocessing import Pool
from datetime import date


def ba_for_haplotype(tup):
    k, keys, i = tup
    max_col = allele_ba_predictions[keys].max(axis=1)
    df = pd.DataFrame(max_col, columns=[k])
    return df


def get_parser():
    """Get parser object for calculating haplotype / genotype ba predictions."""

    parser = argparse.ArgumentParser()
    parser.add_argument('--population', '-pop', dest='population', default='World', type=str,
                        help='Target population. Default "World"')
    parser.add_argument('--peptides', '-pep', dest='peptides', required=True, type=str,
                        help='Preprocessed peptide file with every peptide in a new line.')
    parser.add_argument('--frequencies', '-f', dest='f_data', required=True, type=str,
                        help='(Normalized) haplotype / genotype frequency file.')
    parser.add_argument('--allele_binding_affinities', '-ba', dest='predictions', required=True, type=str,
                        help='Binding affinity file for peptides and single alleles.')
    parser.add_argument('--type', '-t', dest='loci', required=True, type=str,
                        help='Type of data that is processed, e.g. mhc1 for haplotype mhc1; mhc2 for haplotype mhc2;'
                             'mhc1_genotype for mhc1 genotype; mhc2_genotype for mhc2 genotype')
    parser.add_argument('--nworkers', '-n', dest='nworkers', default=1, type=int,
                        help='Number of cores to use in parallel.')

    return parser


def main():
    args = get_parser().parse_args()

    filtered_peptides = read_peptides.main(args.peptides)
    print('Read number of peptides', len(filtered_peptides))

    global allele_ba_predictions
    allele_ba_predictions = pd.read_pickle(args.predictions).loc[filtered_peptides]
    print('Loaded ba predictions', allele_ba_predictions.shape)

    haplotypes = pd.read_pickle(args.f_data).keys()

    if args.loci == 'mhc1':
        loci = ['HLA-A', 'HLA-B', 'HLA-C']
    elif args.loci == 'mhc2':
        loci = ['DRB1', 'HLA-DP', 'HLA-DQ']
    elif args.loci == 'mhc1_genotype':
        loci = ['HLA-A', 'HLA-A', 'HLA-B', 'HLA-B', 'HLA-C', 'HLA-C']
    elif args.loci == 'mhc2_genotype':
        loci = ['DRB1', 'DRB1', 'HLA-DP', 'HLA-DP', 'HLA-DQ', 'HLA-DQ']
    else:
        print('Type unknown! Choose from [mhc1, mhc2, mhc1_genotype, mhc2_genotype].')
        exit(-1)

    print('Start collecting keys')
    hap_list = []
    for i, key in enumerate(haplotypes):
        if i % 10000 == 0:
            print(i, 'of', len(haplotypes))
        key = ['unknown' if x=='nope' else x for x in key]
        keys = list(zip(loci, key))
        if all(k in allele_ba_predictions for k in keys):
            hap_list.append((tuple(key), keys, i))
    print('Processed all keys in genotypes / haplotypes')

    print('Create threads')
    pool = Pool(processes=args.nworkers)
    outputs = pool.map_async(ba_for_haplotype, hap_list)

    pool.close()
    pool.join()

    print('Finished all tasks')
    new_df = pd.concat([out for out in outputs.get()], axis=1)
    print('Finished finding maximum for each genotype / haplotype')

    today = date.today().strftime('%d-%b')
    new_df.to_pickle('../../Gifford_Data/' + today + '_' + args.population + '_predictions_for_' + args.loci +
                     '_filtered_peptides.pkl.gz', compression='gzip')

    return new_df


if __name__ == '__main__':
    main()

    countries = ['Europe',
                 'Northeast Asia',
                 'Oceania', 'East Africa',
                 'East Asia', 'South Africa',
                 'Southeast Asia',
                 'Central Africa', 'North America',
                 'South America', 'West Africa',
                 'West Indies', 'South Asia',
                 'North Africa', 'Southwest Asia']

    # for country in countries:
    #     print(country)
    #     haplotypes = pd.read_pickle('../EvaluateVaccines/data/haplotype_freq_no_ld_' + country + '.pkl')
    #     main(threads, mhc1_ba_predictions, haplotypes, ['HLA-A', 'HLA-B', 'HLA-C'], country, 'mhc1')

    # haplotypes = pd.read_pickle('../EvaluateVaccines/data/haplotype_freq_no_ld_' + countries[8] + '.pkl')
    # main(threads, mhc1_ba_predictions, haplotypes, ['HLA-A', 'HLA-B', 'HLA-C'], countries[8], 'mhc1')
