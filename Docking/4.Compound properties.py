import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import os

# set working directory
os.chdir("./Docking")

# impot data
data = pd.read_csv('Enamine_predictions.csv')

# Calculate molecular weight from SMILES column and add to dataframe
data['Mol'] = data['SMILES'].apply(Chem.MolFromSmiles)
data['MW'] = data['Mol'].apply(Descriptors.ExactMolWt)

# Calculate number of rotatable bonds from SMILES column and add to dataframe
data['NumRotatableBonds'] = data['Mol'].apply(Descriptors.NumRotatableBonds)

# Calculate number of rings from SMILES column and add to dataframe
data['NumRings'] = data['Mol'].apply(Descriptors.RingCount)

#save dataframe to csv
data.to_csv('Enamine_predictions_with_stats.csv', index=False)
