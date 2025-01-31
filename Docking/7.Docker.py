from qsprpred.data.storage.tabular.simple import PandasChemStore
from qsprpred.data.chem.identifiers import InchiIdentifier
from qsprpred.data.chem.standardizers.papyrus import PapyrusStandardizer
from qsprpred.data.storage.tabular.hierarchical import PandasRepresentationStore
from spock.prep import Dimorphite
from spock.storage.tabular import SpockStorage
from spock.storage.tabular import SpockProtein
from spock.docking.vina.cpu_local import VinaDockingCPULocal
import shutil
import pandas as pd
import os

def main():
    # set working directory
    os.chdir("./Docking")

    # List of directories to remove
    directories_to_remove = ['./data/SpockStorage', './data/ChemStorage', './data/ProtomerStorage']

    # remove the directories if they exist
    for directory in directories_to_remove:
        if os.path.exists(directory):
            shutil.rmtree(directory)

    # Create a library of compounds stored as a pandas data frame
    library = PandasChemStore(
        name="ChemStorage",
        path="./data/",
        df=pd.read_csv("subset.csv"),
        standardizer=PapyrusStandardizer(),
        identifier=InchiIdentifier(),
    )
    library.save()

    # make a storage for protomers of molecules in our library
    representation_store = PandasRepresentationStore(
        name="ProtomerStorage",
        path="./data/",
        chem_store=library, # our library from before
    )
    dm = Dimorphite(ph_range=(7,7.5), max_variants=5)
    dm.apply_to_storage(representation_store)
    
    representation_store.save()
    
    # Create Spock storage for ligand poses
    store = SpockStorage(
        name="SpockStorage",
        path="./data/",
        ligand_store=representation_store,
    )

    PROTEIN_NAME = "7vhy-prepped"   # name of the protein
    N_CPUS = os.cpu_count()  # number of cpus to use for docking
    EXHAUSTIVENESS = 8  # Vina exhaustiveness parameter (ideally no lower than 8)
    SEED = 42  # random seed for random operations (the meaning of life)
    PROTEIN_FOLDER = './data/Receptor'

    # creating a docking engine
    docking = VinaDockingCPULocal(
        protein=SpockProtein(
            PROTEIN_NAME,
            props={
                "pdb": open(f'{PROTEIN_FOLDER}/{PROTEIN_NAME}.pdb', 'r').read(),
                "pdbqt": open(f'{PROTEIN_FOLDER}/{PROTEIN_NAME}.pdbqt', 'r').read(),
            }
        ),
        n_cpus=N_CPUS,
        box_spec={
            "center": [-12.2, -16.1, 13.3], # setting the the docking box
            "box_size": [16.5, 20.6, 16.7]
        },
        embed_mols=True,  # set to False if conformers are already generated
        exhaustiveness=EXHAUSTIVENESS,
        seed=SEED,
        timeout=60 # timeout in seconds for each molecule
    )

    # dock the molecules
    store.nJobs = N_CPUS # n_cpus can also be set on the storage itself and the docking engine will use it
    docking.dock_storage(storage=store, chunk_size=1, save=True, overwrite=True)

    # save the storage
    store.save()

    # print the summary of the storage
    # print(store.getSummary())
    # print(representation_store.getSummary())
    # print(representation_store.getSummary()['num_reps'].sum())

if __name__ == "__main__":
    main()