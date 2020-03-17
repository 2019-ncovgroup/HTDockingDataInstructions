import pandas as pd
import pickle
from tqdm import tqdm
import numpy as np

radical = pd.read_pickle("radical.pkl")
balsam = pd.read_pickle("balsam_out.pkl")

print(radical.shape)
print(balsam.shape)

df = pd.merge(left=radical, right=balsam, on='smiles', how='outer', suffixes=('_radical', '_balsam'))

for i in tqdm(range(1, df.shape[1])):
    df.iloc[:, i] = pd.to_numeric(df.iloc[:,i])

cols  = list(df.columns[1:])
cols = list(filter(lambda x : "_radical" in x or "_balsam" in x, cols))
print("overlap on", cols)
cols.sort()
for i in range(0, len(cols), 2):
    tmp = df.loc[:, [cols[i], cols[i+1]]]
    tmp = np.array(tmp)
    df["_".join(cols[i].split("_")[:-1])] = np.nanmean(tmp, axis=1)
    df = df.drop([cols[i], cols[i+1]], axis=1)

print(cols)
print(df.head())
print(df.columns)
print(df.shape)
#
df.to_pickle('out.pkl')
df.to_csv('out.csv', index=False)
df.to_csv('out.csv.gz', index=False)


#prep ML

df.iloc[:, 1:] = np.nan_to_num(np.clip(np.array(df.iloc[:, 1:]), -10000, 0), nan=0)
targets = df.columns[1:].tolist()



keep_targets = []
for target in targets:
    if np.quantile(df.loc[:, target], 0.025) > -8 or np.abs(np.mean(df.loc[:, target]) - np.quantile(df.loc[:, target], 0.025)) < 1:
        print("removing {} for ml".format(target))
    else:
        keep_targets.append(target)

df = df[['smiles'] + keep_targets]
print("Kept:", keep_targets)
print(df.shape)
print(df.head())

df.to_pickle('dock_out_ml_v1.pkl')
df.to_csv('dock_out_ml_v1.csv', index=False)
df.to_csv('dock_out_ml_v1.csv.gz', index=False)
