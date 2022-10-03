import sys


def main(file):
    peptides = []
    for line in open(file, 'r'):
        if line.startswith('#'):
            continue
        pep = line.strip()
        peptides.append(pep)
    return peptides

if __name__ == '__main__':
    f = sys.argv[1]
    main(f)