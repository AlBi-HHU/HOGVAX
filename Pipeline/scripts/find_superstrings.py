import sys
import read_peptides

def main(peptides):
    # sort peptides in reversed order
    sorted_peptides = sorted(peptides, key=len, reverse=True)
    # first read in all peptides and identify superstrings
    superstrings = {}
    for i, peptide in enumerate(sorted_peptides):
        if i % 1000 == 0:
            print(i, 'of', len(sorted_peptides))
        for j in range(i+1, len(sorted_peptides)):
            sec_peptide = sorted_peptides[j]
            if sec_peptide in peptide:
                if peptide in superstrings:
                    superstrings[peptide].append(sec_peptide)
                else:
                    superstrings[peptide] = [sec_peptide]

    return superstrings

if __name__ == '__main__':
    peptide_file = sys.argv[1]
    peptides = read_peptides.main(peptide_file)
    main(peptides)
