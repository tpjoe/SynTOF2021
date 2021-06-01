"""
This script essentially generate xy coordinates for single cell plot for visualization in R
in script 18) using TSNE. (Just because python implementation of TSNE is much faster than R).
"""
from cuml.manifold import UMAP, TSNE
from cuml import KMeans
import cudf
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import random
from glob import glob
import flowkit as fk
import re
from sklearn.manifold import TSNE as skTSNE


# plot clusters ----------------------------------------------------------------------------------------------------
mc = pd.read_csv('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv')
mc.loc[:, 'sample'] = mc.loc[:, 'sample'].apply(lambda x: re.sub('_BC\d+', '', x))

hidden_ln = cudf.read_csv('R_py_exchange/hidden_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv').iloc[:, 1:].to_pandas()
hidden_ln = hidden_ln.loc[hidden_ln.loc[:, 'sample'].apply(lambda x: ('HF14-017.fcs' not in x) & ('HF14-025.fcs' not in x) & ('HF14-083.fcs' not in x)), :]
# hidden_lbd = cudf.read_csv('R_py_exchange/hidden_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1_LBD.csv').iloc[:, 1:].to_pandas()
# hidden_ad = cudf.read_csv('R_py_exchange/hidden_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1_PHAD.csv').iloc[:, 1:].to_pandas()
# hidden_ = pd.concat([hidden_ln, hidden_lbd, hidden_ad], axis=0).reset_index(drop=True)

hidden_ = hidden_ln
regions = hidden_.loc[:, 'sample'].apply(lambda x: x.split('_')[0])
prob = regions.value_counts()/regions.value_counts().max()
p = regions.map({'BA9': 1/0.67088, 'DLCau': 1/0.52534, 'Hipp': 1})

# sample for visual
np.random.seed(0)
ind = np.sort(np.random.choice(hidden_.shape[0], 90000, p=p/sum(p), replace=False))

hidden = hidden_.iloc[ind, :-1] #.to_pandas().sample(n=50000, replace=False).sort_index().iloc[:, :-1]
mc_ = mc.iloc[list(hidden.index), 0].values


umap = skTSNE(n_components=2, perplexity=150, n_jobs=-1, random_state=99, early_exaggeration=120)
umap_xy = umap.fit_transform(hidden)


plot_data = pd.concat([pd.DataFrame(umap_xy), pd.Series(mc_)], axis=1)
plot_data.columns = ['x', 'y', 'c']

ax = plt.figure()
number_of_colors = len(np.unique(mc_))
set.seed(0)
color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]
# sns.set_palette(color)
ax = sns.scatterplot(x='x', y='y', s=2, hue='c', palette=color, linewidth=0, data=plot_data) # palette="Set2"
# ax.get_legend().remove()
plt.savefig('test_cluster_allEvents.png')

# saving for better visualization in R
plot_data.to_csv('R_py_exchange_afterCluster/umap_sc.csv')
plot_data2 = pd.concat([plot_data, pd.Series(ind+1).reset_index(drop=True)], axis=1)
plot_data2.columns = ['x', 'y', 'c', 'ind+1']
plot_data2.to_csv('R_py_exchange_afterCluster/umap_sc_randomLowNo_allRegions.csv')

