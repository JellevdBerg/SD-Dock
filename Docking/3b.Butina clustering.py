import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.ML.Cluster import Butina
from tqdm import tqdm

def read_smiles(file_path):
    """Read SMILES from a CSV file."""
    df = pd.read_csv(file_path)
    return df

def generate_fingerprints(mols):
    """Generate molecular fingerprints."""
    generator = GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(mol) for mol in mols if mol is not None]
    return fps

def cluster_fingerprints(fps, cutoff=0.35):
    """Cluster fingerprints using Butina algorithm."""
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1-x for x in sims])
    clusters = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return clusters

def select_simplest_structure(clusters, mols, df):
    """Select the simplest structure (based on molecular weight) from each cluster."""
    simplest_mols = []
    simplest_ids = []
    for cluster in clusters:
        cluster_mols = [mols[i] for i in cluster if mols[i] is not None]
        cluster_ids = [df.iloc[i]['ID'] for i in cluster if mols[i] is not None]
        simplest_mol = min(cluster_mols, key=lambda mol: Descriptors.MolWt(mol))
        simplest_id = cluster_ids[cluster_mols.index(simplest_mol)]
        simplest_mols.append(simplest_mol)
        simplest_ids.append(simplest_id)
    return simplest_mols, simplest_ids

def main():
    # File path to the CSV file containing SMILES strings
    file_path = './Enamine_predictions.csv'
    
    # Read SMILES from the CSV file
    df = read_smiles(file_path)

    # Convert SMILES to RDKit molecule objects
    mols = [Chem.MolFromSmiles(smiles) for smiles in df['SMILES']]
    
    # Generate fingerprints
    fps = generate_fingerprints(mols)
    
    # Perform Butina clustering
    clusters = cluster_fingerprints(fps, cutoff=0.35)
    
    # Print the number of clusters created
    print(f"Number of clusters created: {len(clusters)}")
    
    # Select the simplest structure from each cluster
    simplest_mols, simplest_ids = select_simplest_structure(clusters, mols, df)
    
    # Convert the simplest molecules back to SMILES
    simplest_smiles = [Chem.MolToSmiles(mol) for mol in simplest_mols]
    
    # Create a DataFrame to store the results
    result_df = pd.DataFrame({'ID': simplest_ids, 'SMILES': simplest_smiles})
    
    # Merge with original DataFrame to get the Probability (active) column
    result_df = result_df.merge(df[['ID', 'Prediction (active)']], on='ID')
    
    # Save the results to a new CSV file
    result_df.to_csv('./simplest_strucures.csv', index=False)
    
    print("Simplest structures have been selected and saved to 'simplest_structures.csv'.")

if __name__ == "__main__":
    main()