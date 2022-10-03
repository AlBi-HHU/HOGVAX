import glob
import re

def get_vax_len(in_file, out_file):
    out = open(out_file, 'w')
    peptide_list = []
    with open(in_file) as file:
        for line in file:
            peptides_, score = line.split('	')
            peptide_list = peptides_.split('_')
            for p in peptide_list:
                out.write(p + '\n')
            break

    concat_seq = ''.join(peptide_list)
    seq_len = len(concat_seq)
    out.write('# Total length ' + str(seq_len))

    print(concat_seq)
    print(seq_len)
    return seq_len


dict_k = {}
for path in glob.glob('/Users/sara/Documents/VaccinesProject/ivp/Gifford_Code/result_mhc2_maxround*'):
    print(path)
    size = re.search(r'maxround(\d+)', path).group(1)
    file = path + '/best_result.txt'
    out_file = path + '/13Sept_optivax_unlinked_peptides.txt'
    dict_k[size] = get_vax_len(file, out_file)

print(dict_k)

with open('../../Pipeline/snakemake_mhc2_call_file.sh', 'w') as f:
    for key in dict_k:
        f.write('snakemake --use-conda --cores 32 --config mhc=\'14Sept_mhcII_' + key + '\' k=' + str(dict_k[key]) + '\n')
