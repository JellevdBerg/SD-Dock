# SD-Dock
For my personal course on docking analyses

Data is stored within the data folder.
All scripts run accordingly if the repo is correctly cloned. no alterations to file paths and directories have to be made for the workflow to function.
Run order of the scripts and notebooks are in order of numbering. 3b (butina clustering) is optional and not required for the workflow unless your machine/server has plentiful memory available.

This repository relies on a very specific set of python packages and versions in order to function properly, the following step by step plan for recreation is advised:
- create conda env with python version 3.12.2
- install vina using "conda install -c conda-forge vina"
- clone the spock repository developed by M. Sicho from "https://github.com/martin-sicho/spock/tree/dev"
- clone the QSPRpred repository developed by the LACDR-DDS group from "https://github.com/CDDLeiden/QSPRpred/tree/feature/representations_storage"
- CD to the repositories and "pip install ." Spock and QSPRpred (Sock has to be installed before QSPRpred for it to function properly)
- run "pip install rdkit==2023.9.1" in the terminal to change the rdkit version to the correct one
- IF the workflow complains about dimorphite: "conda install -c conda-forge dimorphite-dl=1.3.2"

This should result in a workable conda env that can run all the scripts within this repository.