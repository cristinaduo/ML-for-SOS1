# ML-for-SOS1

This repository contains the original codes and data in *Discovery of Novel SOS1 Inhibitors Using Machine Learning*.



## Requirements
* Python 3.7
* Numpy 1.21.6
* Pandas 1.1.3
* RDKit 2019.09.3
* scikit-learn 1.0.2
* matplotlib 3.5.1
* jupyterlab 3.5.3

## Introduction of relevant datasets


### Dataset for model validation
The 'labelled_molecules' folder houses molecules related to SOS1 as documented in the literature. Each molecule's information has been transformed into molecular fingerprints, and these fingerprints are concatenated with the 10% fixed data to construct the external validation set.


## Script Description


#### Data cleaning, model construction, and optimization
`Machine Learning on SOS1 inhibitors.ipynb`

#### External dataset construction for model validation
`Model Validation.ipynb`

#### EGFR-related data virtual screening
`VS of EGFR ChEMBL.ipynb`

#### Utils
`utils.py`, `rdkit_utils.py` and `rdkit_utilies.py`
