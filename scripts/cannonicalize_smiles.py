from rdkit import Chem
from tqdm import tqdm
import pandas as pd

def cannon_smile(smi):
    try:
        cannon = Chem.MolFromSmiles(smi)
        if cannon is not None:
            return Chem.MolToSmiles(cannon, canonical=True)
        else:
            print(smi)
    except:
        print('error', smi)
    return smi

df = pd.read_pickle("dock_out_ml_v1.pkl")
smiles = df.smiles.tolist()
canon_smiles = [cannon_smile(smile) for smile in tqdm(smiles)]
