import pandas as pd
import re

keep = []
indices = []
with open('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/SarsCoVCleavage/NetChop3.1_sars_cov_2_predictions.txt') as file:
    for line in file:
        line = line.lstrip()
        if line.startswith('#') or line.startswith('-'):
            continue
        if line.startswith('pos'):
            indices.append(len(keep))
            line = re.sub(pattern='\s+', repl='\t', string=line)
            line = line.strip().split('\t')
            keep.append(line)
            continue
        if re.match('^\d', line):
            match = re.search(r'\s+(\D+)\n', line)
            if match:
                line = line.replace('Base', '').replace('_', '')
                line = re.sub(pattern='\s+', repl='\t', string=line)
                line = line.strip().split('\t')
                keep.append(line)

print(keep)
df = pd.DataFrame(keep)
df = df.iloc[: , :-1]
print(df)
df.to_pickle('/home/sara/Documents/VaccinesProject/ivp/Proteasomal_Cleavage_Predictions/SarsCoVCleavage/NetChop3.1_sars_cov_2_predictions.pkl')