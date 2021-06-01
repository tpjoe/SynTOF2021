# SynTOF2021
This repository contains code to replicate results from the paper titled "Single-synapse analysis of Alzheimer’s disease implicate pathologic tau, DJ1, CD47, and ApoE".


The descriptions of each scripts and folders are written below:

All scripts are numbered in a sequential order to achieve the presented analyses, these include:
* Script# 1.1-1.3: Bar plots and spider plots of the pseudo-bulk analysis (R and Python)
* Script# 2: Converts the single-synapse data file format from .txt files to .fcs files (R)
* Script# 3: Performs autoencoder clustering of the single synapse data (Python)
* Script# 4: Performs metaclustering of the resulting autoencoder clusters (R)
* Script# 5: Performs post-synaptic reference correction (R)
* Script# 6: Generates the final frequency and expression matrix used for downstream analysis (R)
* Script# 7: Visualizes mean phenotypic expression values of each cluster using heatmap (R)
* Script# 8: Visualizes frequency of each cluster (R)
* Script# 9: Visualizes correlations between features (R)
* Script# 10: Visualizes cluster through GateFinder (R)
* Script# 11: Performs conventional machine learning on the generated expression feature matrix (Python)
* Script# 12: Visualizes feature importance (R)
* Script# 13: Visualizes model prediction values (R)
* Script# 14: Visualizes boxplots of the selected features stratified by diagnosis groups (R)
* Script# 15: Visualizes all significantly different features between AD vs. Control as well as non-resilient AD vs. Control (R)
* Script# 16: Visualizes changes if resilient AD were to be removed (in a candy bar plot) (R)
* Script# 17: Generates x, y coordinates for single cell vissualization using t-SNE, done in Python because faster than R (Python)
* Script# 18: Visualizes single cell data, colored by expression or cluster assignment, using coordinates from Script$ 17 (R)
* Script# 19: Performs similar task to Script# 4, but for GFAP-EAAT- single synapse data (R)
* Script# 20: Performs similar task to Script# 14, but for GFAP-EAAT- single synapse data (R)
* Script# 21: Performs autoencoder clustering with leave-one-out scheme for testing clustering robustness (Python)
* Script# 22: Performs metaclustering of results from Script# 21 and analyzes its similarity to clusters obtained from Script# 4 (R)
