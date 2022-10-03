import sys
import re
import read_peptides, find_superstrings


def get_all_peptides(all_peptides, vax_peptides):
    superstrings_mhc = find_superstrings.main(all_peptides)
    pep_incl_substr = vax_peptides
    for pep in vax_peptides:
        if pep in superstrings_mhc:
            pep_incl_substr = pep_incl_substr + superstrings_mhc[pep]

    print(len(set(pep_incl_substr)))
    return set(pep_incl_substr)


def remove_embedding(peptides, strip):
    stripped = [p[strip:-strip] for p in peptides]
    return stripped
    
    
# use pipeline output containing all chosen peptides
peptide_file = snakemake.input['chosen_peptides']
ilp_peptides = list(set(read_peptides.main(peptide_file)))

# optional for embedding use embedding length as third argument
if snakemake.params['embed_len']:
    print('###################  Get all peptides inc substrings: embedded will be executed! ##################')
    ilp_peptides = remove_embedding(ilp_peptides, snakemake.params['embed_len'])

# use pipeline input of all peptides from which ilp may choose
all_peptides = snakemake.params['pep']

with open(snakemake.output[0], 'w') as file:
    file.write('\n'.join(get_all_peptides(all_peptides, ilp_peptides)))
