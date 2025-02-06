# SD-Dock
For my personal course on docking analyses

Data is stored within the data folder.
All scripts run accordingly if the repo is correctly cloned. no alterations to file paths and directories have to be made for the workflow to function.
Run order of the scripts and notebooks are in order of numbering. 3b (butina clustering) is optional and not required for the workflow unless your machine/server has plentiful memory available.

### Dependencies
All data is stored within the “data” folder. The scripts are designed to run seamlessly as long as the repository is cloned correctly, with no need for modifications to file paths or directories. The execution order of the scripts and notebooks follows a clear, numbered sequence, ensuring the workflow operates smoothly.

> [!IMPORTANT]
>This repository relies on a very specific set of python packages and versions in order to function properly, the following step by step plan for recreation is advised:
1. **Create a virtual environment** with Python 3.12.2
   ```
   create conda -n spock python=3.12.2
   ```
2. **Install Vina** using:  
   ```
   conda install -c conda-forge vina
   ```
3. **Clone the [Spock](https://github.com/martin-sicho/spock/tree/dev) repository** (developed by M. Sicho):  
   ```
   git clone https://github.com/martin-sicho/spock/tree/dev
   ```
4. **Clone the [QSPRpred](https://github.com/CDDLeiden/QSPRpred/tree/feature/representations_storage) repository** (developed by the LACDR-CDD group):  
   ```
   git clone https://github.com/CDDLeiden/QSPRpred/tree/feature/representations_storage
   ```
5. **Navigate to each repository directory** using `cd "INSERT_PATH_TO_REPO"` and install each package with:  
   ```
   pip install .
   ```  
   Make sure **Spock** is installed before **QSPRpred**.
6. **Install the correct version of RDKit**:  
   ```
   pip install rdkit==2023.9.1
   ```
7. **If the workflow reports an issue with Dimorphite**, install the required version:  
   ```
   conda install -c conda-forge dimorphite-dl=1.3.2
   ```

This should result in a workable conda env that can run all the scripts within this repository.