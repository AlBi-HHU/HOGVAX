import sys
import re
import find_superstrings


def get_all_peptides(tree, leaves, vax_peptides, path):
    superstrings_mhc = find_superstrings.aho_corasick_algorithm(tree, leaves, vax_peptides)
    pep_incl_substr = vax_peptides
    for pep in vax_peptides:
        if pep in superstrings_mhc:
            pep_incl_substr = pep_incl_substr + superstrings_mhc[pep]

    print('Number of peptides included in vaccine', len(set(pep_incl_substr)))
    return set(pep_incl_substr)


def remove_embedding(peptides, strip):
    stripped = [p[strip:-strip] for p in peptides]
    return stripped
    

def recreate_unembedded_peptides(input_peptides, hogvax_peptides, embedding_length, path):
    # optional for embedding use embedding length as third argument
    if embedding_length > 0:
        print('###################  Get all peptides inc substrings: embedded will be executed! ##################')
        hogvax_peptides = remove_embedding(hogvax_peptides, embedding_length)

    # use input of all peptides from which ilp may choose
    unembedded_substr_incl_peptides = get_all_peptides(input_peptides, hogvax_peptides)

    with open(path + 'pep_out/' + str(len(input_peptides)) + '_chosen_peptides_hogvax_inc_substrings.txt', 'w') as file:
        file.write('\n'.join(unembedded_substr_incl_peptides))

    with open(path + 'pep_out/' + str(len(input_peptides)) + '_hogvaxine.txt', 'a') as out:
        out.write('> MHC optimized combined peptide vaccine sequence concatenated incl all substring peptides\n')
        out.write(''.join(unembedded_substr_incl_peptides))
