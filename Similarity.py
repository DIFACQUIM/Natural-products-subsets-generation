
from tkinter import E
import pandas as pd
import numpy as np
from rdkit.Chem import MACCSkeys
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import pdist

# Read files
subset_a = pd.read_csv("/home/ana/Documentos/FINGERPRINT/HANNA/GitHub/MaxMin/DATASET/UNPD_subset_A_14994_curada_NPL_score.csv")
subset_b = pd.read_csv("/home/ana/Documentos/FINGERPRINT/HANNA/GitHub/MaxMin/DATASET/UNPD_subset_B_7497_curada_NPL_score.csv") 
subset_c = pd.read_csv("/home/ana/Documentos/FINGERPRINT/HANNA/GitHub/MaxMin/DATASET/UNPD_subset_C_4998_curada_NPL_score.csv") 
subset_real = pd.read_csv("/home/ana/Documentos/FINGERPRINT/HANNA/GitHub/MaxMin/DATASET/REAL_subset_3809_curada_NPL_score.csv")
biofacquim = pd.read_csv("/home/ana/Documentos/FINGERPRINT/HANNA/GitHub/MaxMin/DATASET/BIOFACQUIM.V2_NPL_Score.csv")
dnmt1 = pd.read_csv("/home/ana/Documentos/FINGERPRINT/HANNA/GitHub/MaxMin/DATASET/DNMT1_ACTIVITIES_714_27SEP22_NPL_score.csv")

smi_a = list(subset_a["SMILES_chiral"])
smi_b = list(subset_b["SMILES_chiral"])
smi_c = list(subset_c["SMILES_chiral"])
smi_real = list(subset_real["SMILES_chiral"])
smi_biofacquim = list(biofacquim["SMILES_chiral"])
smi_dnmt1 = list(dnmt1["SMILES_chiral"])

print(len(smi_biofacquim))
print(len(smi_dnmt1))

def ECFP (smi, r):
    fps = pd.DataFrame([[int(y) for y in AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), r, nBits=1024).ToBitString()] for x in smi])
    SimMat = 1 - pdist(fps[[x for x in range(1024)]], metric="jaccard") # Similarity Matrix
    #print(SimMat.shape)
    SimMat = round(np.median(SimMat), 3) 
    return SimMat
def MACCS (smi):
    fps = pd.DataFrame([[int(y) for y in MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(x)).ToBitString()] for x in smi])
    SimMat = 1 - pdist(fps[[x for x in range(167)]], metric="jaccard") # Similarity Matrix
    print(SimMat.shape)
    SimMat = round(np.median(SimMat), 3) 
    return SimMat

print(ECFP(smi_biofacquim, 2))
print(ECFP(smi_dnmt1, 2))
print(MACCS(smi_dnmt1))

ecfp_2 = [ECFP(smi_a, 2), ECFP(smi_b, 2), ECFP(smi_c, 2), ECFP(smi_real, 2), ECFP(smi_biofacquim, 2), ECFP(smi_dnmt1, 2)]
print(ecfp_2)
ecfp_3 = [ECFP(smi_a, 3), ECFP(smi_b, 3), ECFP(smi_c, 3), ECFP(smi_real, 3), ECFP(smi_biofacquim, 3), ECFP(smi_dnmt1, 3)]
print(ecfp_3)
MACCS_keys = [MACCS(smi_a), MACCS(smi_b), MACCS(smi_c), MACCS(smi_real), MACCS(smi_biofacquim), MACCS(smi_dnmt1)]
print(ecfp_3)
Collection = ["SUBSET_A", "SUBSET_B", "SUBSET_C", "SUBSET_REAL", "BIOFACQUIM", "DNMT1"]
Compounds = [len(subset_a), len(subset_b), len(subset_c), len(subset_real), len(biofacquim), len(dnmt1)]

arr = np.array([Collection, Compounds, MACCS_keys, ecfp_2, ecfp_3])
arr = np.transpose(arr)
print("array:", arr)

# Dataframe
columns = ["DATA SET", "Compounds", "MACCS_keys(median)", "ECFP_4 (median)", "ECFP_6 (median)"]
fps = pd.DataFrame(arr, columns=columns)
print(fps)

# Save Dataframe
fps.to_csv("Fingerprints_median_similarity.csv", sep=",", index=False)