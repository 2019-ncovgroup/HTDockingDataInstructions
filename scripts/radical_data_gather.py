import glob
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm


dfs = []
files = glob.glob(sys.argv[1] + "*.out")
for file in tqdm(files):
    recept_name = "_".join(file.split("/")[-1].split("_")[:-1])
    df = pd.read_csv(file, header=None)
    df = df.iloc[:, [2,3]]
    df.columns= ['smiles', 'dock']
    df['receptor'] = recept_name
    dfs.append(df)

df = pd.concat(dfs)
df = df[['smiles', 'receptor', 'dock']]
df = df.pivot_table(values='dock', index='smiles', columns='receptor', aggfunc='first')
df = df.reset_index()
cols = df.columns.tolist()
cols = ['smiles']  + [s + "_dock" for s in cols[1:]]
df.columns = cols

print(df.head())

print(df.shape)

df.to_pickle("radical.pkl")