import argparse
import os
import gzip
import pickle
import networkx as nx
import pandas as pd
import gurobipy as gp
import aho_corasick_trie
import linear_time_hog
import find_superstrings
import filter_substrings_in_ba_predictions
import ba_for_embedded_peptides
import hogvax_ilp
import draw_hog
import get_all_peptides_incl_substrings


def get_parser():
    """Get parser object for combinatorial vaccine design."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--k', '-k', required=True, dest='k', type=int, help='Maximal length of vaccine sequence')
    parser.add_argument('--populations', '-pop', dest='populations', default=['World'], nargs='+',
                        help='Target population(s). Default "World"')
    parser.add_argument('--peptides', '-pep', dest='peptides', required=True, type=str,
                        help='Preprocessed peptide file with every peptide in a new line.')
    parser.add_argument('--allele-frequencies', '-af', dest='f_data', required=True, type=str,
                        help='(Normalized) allele frequency file.')
    parser.add_argument('--ba-threshold', '-t', dest='ba_threshold', default=0.638, type=float,
                        help='If provided, binding affinities are converted to binary data.')
    parser.add_argument('--binding-affinities', '-ba', dest='ba_matrix', required=True, type=str,
                        help='Binding affinity file for input peptides and alleles.')
    parser.add_argument('--required_epitopes', '-epi', dest='required_epitopes',
                        help='File of peptides you want to be present in vaccine')
    parser.add_argument('--min-hits', '-mh', dest='min_hits', default=1, type=int,
                        help='Minimum number of hits for an allele to be covered')
    parser.add_argument('--maximize-peptides', dest='maximize_peptides',  action='store_const', const=True,
                        default=False, help='Maximize number of peptides in the vaccine in a second optimization')
    parser.add_argument('--embedding-length', default=0, type=int, help='Set length of embedding if used')
    parser.add_argument('--embedded-peptides', type=str, help='File containing embedded peptides')
    parser.add_argument('--embedded-epitope_features', type=str, help='Path to embedded epitope features')
    parser.add_argument('--outdir', '-o', default='', type=str, help='Output directory')
    parser.add_argument('--verbose', '-v', dest='logging_enabled', action='store_true')

    return parser


def binarize_entries(df, th):
    df[df >= th] = 1
    df[df < th] = 0
    df = df.astype(int)
    return df


def read_peptides(file):
    peptides = []
    for line in open(file, 'r'):
        if line.startswith('#'):
            continue
        pep = line.strip()
        peptides.append(pep)
    return peptides


def log(string):
    if logging_enabled:
        print(string)


def main():
    args = get_parser().parse_args()

    global logging_enabled
    logging_enabled = args.logging_enabled

    if args.outdir:
        if not args.outdir.endswith('/'):
            args.outdir = args.outdir + '/'
        if not os.path.exists(args.outdir):
            os.mkdir(args.outdir)

    # normal peptides
    log('Reading peptides')
    peptides = read_peptides(args.peptides)
    pep_count = len(peptides)

    required_peptides = []
    if args.required_epitopes:
        log('Reading required epitopes')
        required_peptides = read_peptides(args.required_epitopes)

    drawing_enabled = False
    if pep_count < 30:
        log('Number of peptides below 30 -> drawing enabled')
        drawing_enabled = True

    log('Read frequencies')
    f_data = pd.read_pickle(args.f_data)
    log('Read ba predictions')
    ba_matrix = pd.read_pickle(gzip.open(args.ba_matrix))
    log('Binarize ba predictions')
    bin_matrix = binarize_entries(ba_matrix, args.ba_threshold)

    log('Build Aho-Corasick Trie')
    # needed for identification of superstrings
    ac_leaves_dict, ac_trie = aho_corasick_trie.main(peptides, args.outdir, logging_enabled, drawing_enabled)
    log('Identify substrings')
    # first read in all peptides and identify superstrings
    superstrings = find_superstrings.aho_corasick_algorithm(ac_trie, ac_leaves_dict, peptides)
    log('Modify ba predictions for substrings')
    # use unembedded peptides for substr modifications
    mod_ba_df = filter_substrings_in_ba_predictions.modify_ba_predictions(superstrings, bin_matrix, pep_count, args.outdir)
    # make pandas df sparse
    mod_ba_df = mod_ba_df.astype(pd.SparseDtype(int))

    # if using embedded peptides
    if args.embedding_length > 0:
        log('Modify ba predictions for embedding')
        # use substr modified ba data and reset index to embedded peptides
        embedded_peptides = read_peptides(args.embedded_peptides)
        embed_count = len(embedded_peptides)
        mod_ba_df = ba_for_embedded_peptides.embedding_ba_predictions(args.embedded_epitope_features,
                                                                      mod_ba_df,
                                                                      embedded_peptides,
                                                                      args.outdir)
        log('Build Aho-Corasick Trie for embedded peptides')
        # needed for identification of superstrings
        leaves_dict, ac_trie = aho_corasick_trie.main(embedded_peptides, args.outdir, logging_enabled, drawing_enabled)
        log('Create HOG with embedded peptides')
        leaves_dict, hog = linear_time_hog.compute_hog(str(embed_count), embedded_peptides, args.outdir, logging_enabled, drawing_enabled)
    else:
        log('Create HOG')
        leaves_dict, hog = linear_time_hog.compute_hog(str(pep_count), ac_leaves_dict, ac_trie.copy(), args.outdir,
                                                       logging_enabled, drawing_enabled)

    full_known_alleles = gp.tuplelist(key for key in f_data.keys() if key in mod_ba_df.keys())

    log('Start HOG ILP')
    hogvax_peptides = hogvax_ilp.hogvax(k=args.k,
                                        alleles=full_known_alleles,
                                        freq_vector=f_data,
                                        B_matrix=mod_ba_df,
                                        leaves=leaves_dict,
                                        pep_count=str(pep_count),
                                        graph=hog,
                                        min_hits=args.min_hits,
                                        populations=args.populations,
                                        path=args.outdir,
                                        optional_peptides=required_peptides,
                                        maximize_peptides=args.maximize_peptides,
                                        logging=args.logging_enabled,
                                        coloring=drawing_enabled)

    log('Get all chosen peptides and their substrings')
    # recreate (unembedded) input peptides and get all their substrings
    get_all_peptides_incl_substrings.get_all_peptides(ac_trie, ac_leaves_dict, hogvax_peptides, args.outdir)


if __name__ == '__main__':
    main()
