# NetGPT
**NetGPT**: **N**etwork-based **G**enotype-**P**henotype **T**ransformation identifies key signaling crosstalks for prognosis in pan-cancer study

![Figure1](https://github.com/user-attachments/assets/dd9163be-850f-4732-b5ae-40c6eef7a2f4)


![figure4A](https://github.com/4to1stfloor/NetGPT/assets/115065099/41d94253-675d-497f-9bdb-ce2033446f18)


## Repository Structure

This repository contains the following files and directories:

- `NetGPT/`: Directory containing source code
  - `model.R`: Script for data uploading, data preprocessing, training the model, hyperparameter tunning, and save the results
  - Will be added 
- `data/`: Directory containing input datasets
  - `TCGA-XXXX_pathwayeach_all_log.rds`: merged tcga transcriptome data and somatic mutation data as pathway levels
  - `TCGA-XXXX_pathwayeach_all_log.rds`: merged tcga transcriptome data and somatic mutation data as pathway-link levels
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

## Usage

Instructions on how to use the code or reproduce the results.

## To-do-list
- [x] Upload the pretrained model
- [ ] Provide the pretraining code
- [ ] Examples for multi-omics integration, networkpropagation, perturbation prediction
- [ ] Example code for predicting prognosis
- [ ] Refactoring code
