import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import DataStructs
from sys import argv

data = pd.read_csv(argv[1], sep=",")
print(len(data))
# smi selected # n=1
smi_selected = data.sample(n=1, random_state=42)
smi_selected.reset_index(drop=True, inplace=True)
print(smi_selected.shape)

# smi from selected # n = len(data)-1
smi_from_selected = data.drop(smi_selected.index)
smi_from_selected.reset_index(drop=True, inplace=True)
print(smi_from_selected.shape)


def MaxMin(smi_selected, smi_from_selected):
    print(len(smi_selected))
    while len(smi_selected)<500: # Number of compounds in subset
        fps_0 = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x),
                                             2, nBits=1024) for x in smi_selected]
        fps_1 = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x),
                                             2, nBits=1024) for x in smi_from_selected]
        qu, ta, sim, pair = [], [], [], []

        #Compare all fingerprint (fp) pairwise without duplicates.
        for i,j in zip(range(len(fps_0)), range(len(fps_1))):
            s = DataStructs.BulkTanimotoSimilarity(fps_0[i], fps_1[j:])
        #Collect the SMILES and values
        #Each  value (m) from Tanimoto similarity 
            for m in range(len(s)):
                qu.append(smi_selected[i])
                ta.append(smi_from_selected[j:][m])
                sim.append(s[m]) # lista
        # lowest similarity value
        similitud= sim
        min_value = None

        for value in similitud:
            if (min_value is None or value < min_value):
                min_value=value
                min_value_index = similitud.index(min_value)
        # Delate SMILES selecting
        smi_from_selected.remove(str(ta[min_value_index])) # Eliminar smiles seleccionado
        #Add new SMILES to subset
        smi_selected.append(ta[min_value_index])
        print("le", len(smi_selected))
    return smi_selected, smi_from_selected, qu, ta, sim

SUB_SET = MaxMin(list(smi_selected["SMILES_chiral"]), list(smi_from_selected["SMILES_chiral"]))
SUB_SET = SUB_SET[0]
SUB_SET = pd.DataFrame(data=SUB_SET, columns=["SMILES_chiral"])
SUB_SET.to_csv(argv[1][:-4] + "MaxMin.csv", sep=",", index=False)

