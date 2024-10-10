# NetG2P
[![R](https://img.shields.io/badge/R-%3E%3D%204.3.1-blue)](https://www.r-project.org/)

**NetG2P**: **N**etwork-based **G**enotype-to**2**-**P**henotype transformation identifies key signaling crosstalks for prognosis in pan-cancer study

## Graphical abstract
![Figure1](https://github.com/user-attachments/assets/dd9163be-850f-4732-b5ae-40c6eef7a2f4)


![figure4A](https://github.com/4to1stfloor/NetGPT/assets/115065099/41d94253-675d-497f-9bdb-ce2033446f18)

## Usage

For producing network propagation with GC clustering method, follow the code below.

```bash
$ git clone this-repo-url
$ cd ./NetGPT/NetGPT/
$ Rscript Netgpt_GPTmodule.R --expression your expression file --mutation your mutation file --output your path
```

## Repository Structure

This repository contains the following files and directories:

- `NetGPT/`: Directory containing source code
  - `model.R`: Script for data uploading, data preprocessing, training the model, hyperparameter tunning, and save the results
  - Will be added 
- `data/`: Directory containing input datasets
  - `TCGA-XXXX_pathwayeach.rds`: merged tcga transcriptome data and somatic mutation data as pathway levels
  - `TCGA-XXXX_pathwaylink.rds`: merged tcga transcriptome data and somatic mutation data as pathway-link levels
  - `TCGA-XXXX_dual_with_duration.rds`: merged `TCGA-XXXX_pathwayeach.rds`, `TCGA-XXXX_pathwaylink.rds` and duration data from TCGA clinical data
- `tests/`: Will be added 
- `tutorials/`: Will be added 

## File Descriptions

### `NetGPT/model.R`
This script train the model by preprocessed TCGA data. It performs the following steps:
1. Loads preprocessed TCGA data
2. Cleans and split the data.
3. Feature selection
4. Model training
5. Performance evaluation
6. Retraining 3-5 steps until are met the criteria while ablating the noise features
7. Outputs (best model and metrics) to the `your path` directory

### `NetGPT/analysis.R`
will be added. It includes:
- will be added.

## To-do-list
- [x] Upload the pretrained model
- [x] Upload the preprocessed input data
- [x] Upload the input data with clinical information
- [x] Provide the pretraining code for GPT module
- [ ] Provide the pretraining code for COF module
- [ ] Examples for multi-omics integration, networkpropagation, perturbation prediction
- [ ] Example code for predicting prognosis
- [ ] Refactoring code
