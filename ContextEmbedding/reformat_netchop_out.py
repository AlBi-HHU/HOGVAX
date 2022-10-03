import re
import os

prefix = '../../Proteasomal_Cleavage_Predictions/24SeptEmbeddingPredictions/MHCI/'
files = [f for f in os.listdir(prefix) if re.search(r'NetChop3\.1_17Aug_mhcI+(_\d+embedded)*_hogvax(_with_mutations)*\.txt$', f)]
# files = ['NetChop3.1_17Aug_mhcI_1embedded_hogvax.txt', 'NetChop3.1_17Aug_mhcI_3embedded_hogvax.txt',
#          'NetChop3.1_17Aug_mhcI_4embedded_hogvax.txt', 'NetChop3.1_17Aug_mhcI_5embedded_hogvax.txt',
#          'NetChop3.1_17Aug_mhcI_8embedded_hogvax.txt']
# prefix = '../../Proteasomal_Cleavage_Predictions/SarsCoVCleavage/'
print(files)

for file in files:
    keep = []
    indices = []
    with open(prefix + file, 'r') as f:
        for line in f:
            line = line.lstrip()
            if line.startswith('#') or line.startswith('-'):
                continue
            if line.startswith('pos'):
                indices.append(len(keep))
                line = re.sub(pattern='\s+', repl='\t', string=line)
                line = line.strip() + '\n'
                keep.append(line)
                continue
            if re.match('^\d', line):
                match = re.search(r'\s+(\D+\d*\D*)(\s*_*\t*)\n', line)
                if match:
                    ident = match.group(1)
                    line = line.replace(match.group(1), '"' + ident + '"')
                    line = re.sub(pattern='\s+', repl='\t', string=line)
                    line = line.strip() + '\n'
                    match2 = re.search(r'(\t_\t\D*)\"', line)
                    if match2:
                        line = line.replace(match2.group(1), '')
                    keep.append(line)

    if len(indices) == 2:
        with open('../../Proteasomal_Cleavage_Predictions/24SeptEmbeddingPredictions/MHCI/' + file.replace('.txt', '_reformat.txt'), 'w') as out:
            out.write(''.join(keep[indices[0]:indices[1]]))

        with open('../../Proteasomal_Cleavage_Predictions/24SeptEmbeddingPredictions/MHCI/' + file.replace('.txt', '_reformat_concat.txt'), 'w') as out:
            out.write(''.join(keep[indices[1]:-1]))
    else:
        with open('../../Proteasomal_Cleavage_Predictions/SarsCoVCleavage/' + file.replace('.txt', '_reformat.txt'), 'w') as out:
            out.write(''.join(keep[0]))
            for i in range(len(indices)):
                end = i+1
                if end >= len(indices):
                    end = -1
                out.write(''.join(keep[indices[i]+1:indices[end]]))
