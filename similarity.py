
from tkinter import E
import pandas as pd
import numpy as np
from rdkit.Chem import MACCSkeys
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import pdist

# Read files
smi_a = pd.read_csv("UNPD_subset_A_14994_curada_NPL_score.csv")
smi_b = pd.read_csv("UNPD_subset_B_7497_curada_NPL_score.csv") 
smi_c = pd.read_csv("UNPD_subset_C_4998_curada_NPL_score.csv") 
smi_unpd = pd.read_csv("UNPD_NPL_Score.csv")
smi_biofacquim = pd.read_csv("BIOFACQUIM.V2_NPL_Score.csv")
smi_dnmt1 = pd.read_csv("DNMT1_ACTIVITIES_714_27SEP22_NPL_score.csv")

smi_a = list(smi_a["SMILES_chiral"])
smi_b = list(smi_b["SMILES_chiral"])
smi_c = list(smi_c["SMILES_chiral"])
unpd = smi_unpd[["SMILES_chiral"]]
smi_biofacquim = list(smi_biofacquim["SMILES_chiral"])
smi_dnmt1 = list(smi_dnmt1["SMILES_chiral"])

print(len(smi_biofacquim))
print(len(smi_dnmt1))

# sample 1
print(unpd.shape)
unpd1 = unpd.sample(n=1000, random_state=42)
unpd1.reset_index(drop=True, inplace=True) # Reset index
print(unpd1.shape)
unpd = unpd.drop(unpd1.index) # Remove unpd1 from unpd
unpd.reset_index(drop=True, inplace=True) # Reset index
print(unpd.shape)
# sample2
unpd2 = unpd.sample(n=1000, random_state=42)
unpd2.reset_index(drop=True, inplace=True) # Reset index
print(unpd2.shape)
unpd = unpd.drop(unpd2.index) # Remove unpd1 from unpd
unpd.reset_index(drop=True, inplace=True) # Reset index
print(unpd.shape)
# sample3
unpd3 = unpd.sample(n=1000, random_state=42)
unpd3.reset_index(drop=True, inplace=True) # Reset index
print(unpd3.shape)
unpd = unpd.drop(unpd3.index) # Remove unpd1 from unpd
unpd.reset_index(drop=True, inplace=True) # Reset index
print(unpd.shape)
# sample4
unpd4 = unpd.sample(n=1000, random_state=42)
unpd4.reset_index(drop=True, inplace=True) # Reset index
print(unpd4.shape)
unpd = unpd.drop(unpd4.index) # Remove unpd1 from unpd
unpd.reset_index(drop=True, inplace=True) # Reset index
print(unpd.shape)
# sample5
unpd5 = unpd.sample(n=1000, random_state=42)
unpd5.reset_index(drop=True, inplace=True) # Reset index
print(unpd5.shape)
unpd = unpd.drop(unpd5.index) # Remove unpd1 from unpd
unpd.reset_index(drop=True, inplace=True) # Reset index
print(unpd.shape)

unpd_sample= pd.concat([unpd1, unpd2, unpd3, unpd4, unpd5])
print(unpd_sample.shape)
unpd_sample = unpd_sample.reset_index(drop=True)
smi_unpd = list(unpd_sample["SMILES_chiral"])

def ECFP (smi, r):
    fps = pd.DataFrame([[int(y) for y in AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), r, nBits=1024).ToBitString()] for x in smi])
    SimMat = 1 - pdist(fps[[x for x in range(1024)]], metric="jaccard") # Similarity Matrix
    #print(SimMat.shape)
    SimMat = round(np.median(SimMat), 3) 
    return SimMat
def MACCS (smi):
    fps = pd.DataFrame([[int(y) for y in MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(x)).ToBitString()] for x in smi])
    SimMat = 1 - pdist(fps[[x for x in range(167)]], metric="jaccard") # Similarity Matrix
    #print(SimMat.shape)
    SimMat = round(np.median(SimMat), 3) 
    return SimMat

def ECFP_UNPD_sample(dataframe, radious):
    for i in range(10):
        fp = []
        sample = dataframe.sample(2000, random_state=i).copy()
        smi = list(sample["SMILES_chiral"])
        fps = ECFP(smi, radious)
        fp.append(fps)
    return round(np.median(fp), 3)

def MACCSkeys_UNPD_sample(dataframe):
  for i in range(10):
    fp = []
    sample = dataframe.sample(2000, random_state=1).copy()
    smi = list(sample["SMILES_chiral"])
    fps = MACCS(smi)
    fp.append(fps)
  return round(np.median(fp), 3)

print(ECFP(smi_biofacquim, 2))
print(ECFP(smi_dnmt1, 2))
print(MACCS(smi_dnmt1))

ecfp_2 = [ECFP(smi_a, 2), ECFP(smi_b, 2), ECFP(smi_c, 2), ECFP(smi_unpd, 2), ECFP(smi_biofacquim, 2), ECFP(smi_dnmt1, 2)]
print(ecfp_2)
ecfp_3 = [ECFP(smi_a, 3), ECFP(smi_b, 3), ECFP(smi_c, 3), ECFP(smi_unpd, 3), ECFP(smi_biofacquim, 3), ECFP(smi_dnmt1, 3)]
print(ecfp_3)
MACCS_keys = [MACCS(smi_a), MACCS(smi_b), MACCS(smi_c), MACCS(smi_unpd), MACCS(smi_biofacquim), MACCS(smi_dnmt1)]
print(ecfp_3)
Collection = ["UNPD_SUBSET_A", "UNPD_SUBSET_B", "UNPD_SUBSET_C", "UNPD", "BIOFACQUIM", "DNMT1"]
Compounds = [len(smi_a), len(smi_b), len(smi_c), len(smi_unpd), len(smi_biofacquim), len(smi_dnmt1)]

arr = np.array([Collection, Compounds, MACCS_keys, ecfp_2, ecfp_3])
arr = np.transpose(arr)
print("array:", arr)

# Dataframe
columns = ["DATA SET", "Compounds", "MACCS_keys(median)", "ECFP_4 (median)", "ECFP_6 (median)"]
fps = pd.DataFrame(arr, columns=columns)
print(fps)

# Save Dataframe
fps.to_csv("Fingerprints_median_similarity.csv", sep=",", index=False)
