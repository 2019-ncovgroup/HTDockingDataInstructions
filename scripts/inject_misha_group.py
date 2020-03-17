import pickle
import numpy as np
import pandas as pd
import sys

def parse_ver(s):
    return s.split("/")[-2]

def parse_recept(s):
    fname = "_".join(s.split("/")[-1].split("_")[:-1])
    return fname

with open(sys.argv[1], 'rb') as f:
    data = pickle.load(f)

df = pd.DataFrame(data)
df.columns = ['dock', 'smiles', 'smiles_dbase', 'receptor_file']

df['receptor'] = df.iloc[:, -1].apply(parse_recept)

# filter out old things not needed
df = df[df.iloc[:, -2].apply(parse_ver) == 'receptorsV2']

#convert nan to zeros
df.dock = np.nan_to_num(df.dock, nan=0)

df = df[['smiles', 'receptor', 'dock']]

df = df.pivot_table(values='dock', index='smiles', columns='receptor', aggfunc='first')

df = df.reset_index()
cols = df.columns.tolist()
cols = ['smiles']  + [s + "_dock" for s in cols[1:]]
df.columns = cols
print(df.head())
df.to_pickle("balsam_out.pkl")
# df.to_csv("misha_out.csv")