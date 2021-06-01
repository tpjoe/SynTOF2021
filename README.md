# SynTOF2021
## This repository contains code to replicate results from the paper titled "Single-synapse analysis of Alzheimer’s disease implicate pathologic tau, DJ1, CD47, and ApoE".  <br>
Synaptic molecular characterization is limited for Alzheimer’s disease (AD). We used mass cytometry to quantify 38 probes in approximately 17 million single synaptic events from human brains without pathologic change or with pure AD or Lewy body disease (LBD), non-human primates (NHP), and PS/APP mice. Synaptic molecular integrity in humans and NHP was similar. Although not detected in human synapses, Aβ was in PS/APP mice synaptic events. Clustering and pattern identification of human synapses showed expected disease-specific differences, like increased hippocampal pathologic tau in AD and reduced caudate dopamine transporter in LBD, and revealed novel findings including increased hippocampal CD47 and lowered DJ1 in AD and higher ApoE in AD dementia. Our results were supported by multiplex ion beam imaging of intact tissue.

> Currently our manuscript is submitted to journal and is not published yet.

The following codes allow replication of the results, but note that they require the original single synapse data, which currently can be obtained through contacting us at tpjoe_at_stanford.edu. We are working on uploading the data to a public repository.

## The description of each folder is written below:  <br>
* Folder "rss_plots": Contains rss plots for optimal number of clusters during autoencoder clustering
* Folder "R_py_exchange": Contains preprocessed data, mainly autoencoder clustering results and its hidden representation, in a csv format, so that it can be used by R for metaclustering and visualizations
* Folder "df_pre_save": Contains the single synapse data after signal correction by subtracting signals that are known to be non-existent in post-synaptic events.
* Folder "R_py_exchange_afterCluster": Contains preprocessed data, mainly conventional machine model results from python, in a csv format, so that it can be used by R for visualizations
* Folder "figures": Contains generated figures

## The description of the scripts: <br>
Note: All scripts are numbered in a sequential order to achieve the presented analyses

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
* Script# 19: Performs similar tasks to Script# 4, but for GFAP-EAAT- single synapse data (R)
* Script# 20: Performs similar tasks to Script# 14, but for GFAP-EAAT- single synapse data (R)
* Script# 21: Performs autoencoder clustering with leave-one-out scheme for testing clustering robustness (Python)
* Script# 22: Performs metaclustering of results from Script# 21 and analyzes its similarity to clusters obtained from Script# 4 (R)
