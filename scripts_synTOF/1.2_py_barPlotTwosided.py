"""
This script plot the two-sided bar plot in figure one for comparing PHF-tau high vs low characteristics for pseudo-bulk analysis.
"""

# import lib
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('classic')
import numpy as np
import re
import seaborn as sns
# pd.options.display.max_rows = 999
import matplotlib


###################
# For Human data
###################

# import data ------------------------------------------------------------
df = pd.read_csv('../raw_data/Toms_mean_intensity_change2.csv')

# change data type a bit
for i in range(1, 10):
    df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)] = df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)].apply(lambda x: re.sub(',', '', x))
    df.loc[:, 'Sample'+str(i)] = df.loc[:, 'Sample'+str(i)].astype('float')

# do human pre
df_human_pre = df.loc[(df.HuMu=='Hu') & (df.PrePost=='Pre'), :]
df_human_pre = pd.wide_to_long(df_human_pre, ['Sample'], i=['Marker', 'PrePost', 'Region', 'Group', 'HuMu', 'Compare'], \
    j='sample').reset_index()

# Define the sorter
sorter = ['Tau', 'PHF-tau', 'CD56', 'SNAP25', 'CD47', 'VGLUT', 'GAD65', 'GATM']
# Create the dictionary that defines the order for sorting
sorterIndex = dict(zip(sorter, range(len(sorter))))
df_human_pre['Marker_rank'] = df_human_pre['Marker'].map(sorterIndex)
XX = df_human_pre.dropna().sort_values(['Marker_rank', 'Group', 'Marker'], ascending=[True, False, True])

# start plotting
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(5,4))
ax2.yaxis.set_tick_params(labelright=True)
ax1.yaxis.set_tick_params(labelleft=False)
ax1.spines['right'].set_visible(False)
ax1.tick_params(axis='y', which='both', bottom=False)
ax2.spines['left'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(axis='y', which='both', left=False, right=False)       
ax2.tick_params(axis='y', which='both', left=False, right=False)   
ax2.tick_params(axis='x', which='both', top=False) 
ax1.tick_params(axis='x', which='both', top=False)    


ax2.set_xlim(0,60)
ax2.set_xticks(np.arange(0, 61, 20))
ax1.set_xlim(14000, 21000)
ax1.set_xticks(np.arange(14000, 21001, 2000))
# bars1 = XX.plot(ax=ax1, kind='barh')
# bars2 = XX.plot(ax=ax2, kind='barh')


sns.barplot(ax=ax2,
    data=XX, x="Sample", y="Marker", hue="Group", ci='se', 
    palette=['black', '#fa164f'], capsize=.2)

sns.barplot(ax=ax1,
    data=XX, x="Sample", y="Marker", hue="Group", ci='se', 
    palette=['black', '#fa164f'], capsize=.2)


ax2.legend_.remove()
ax1.legend(loc='lower left').set_title('')
ax1.set_ylabel('')
ax2.set_ylabel('')
ax1.set_xlabel('') 

for tick in ax1.get_xticklabels():
    tick.set_rotation(90)

for tick in ax2.get_xticklabels():
    tick.set_rotation(90)

for tick in ax2.get_yticklabels():
    tick.set_ha('center')

# apply offset transform to all x ticklabels.
dx = 0.45; dy = 0. 
offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
for label in ax2.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)


d = .03
kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
ax2.plot((-d, +d), (-d, +d), **kwargs)      
# ax1.plot((-d, +d),(1 - d, 1 + d), **kwargs)

kwargs.update(transform=ax1.transAxes)  
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  
# ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

ax1.invert_xaxis()
ax2.invert_xaxis()

# plt.tight_layout()
plt.subplots_adjust(right=0.9-0.1, bottom=0.1+0.12, left=0.125-0.05)
plt.rc('font', size=15)
plt.rc('legend', fontsize=15) 
plt.rc('xtick', labelsize=15)
ax2.set_xlabel('% Change PHF-tau$\mathregular{^{high}}$\nvs. PHF-tau$\mathregular{^{low}}$ in Presynaptic', fontsize=13) 
plt.savefig('figures/barplots/hu_pre.pdf')



# for right hand side
# import data ------------------------------------------------------------
df = pd.read_csv('../raw_data/Toms_mean_intensity_change2.csv')

# change data type a bit
for i in range(1, 10):
    df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)] = df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)].apply(lambda x: re.sub(',', '', x))
    df.loc[:, 'Sample'+str(i)] = df.loc[:, 'Sample'+str(i)].astype('float')

# do human pre
df_human_pre = df.loc[(df.HuMu=='Hu') & (df.PrePost=='Post'), :]
df_human_pre = pd.wide_to_long(df_human_pre, ['Sample'], i=['Marker', 'PrePost', 'Region', 'Group', 'HuMu', 'Compare'], \
    j='sample').reset_index()


# Define the sorter
sorter = ['Tau', 'PHF-tau', 'CD56', 'SNAP25', 'CD47', 'VGLUT', 'GAD65', 'GATM']
# Create the dictionary that defines the order for sorting
sorterIndex = dict(zip(sorter, range(len(sorter))))
df_human_pre['Marker_rank'] = df_human_pre['Marker'].map(sorterIndex)
XX = df_human_pre.dropna().sort_values(['Marker_rank', 'Group', 'Marker'], ascending=[True, False, True])

# start plotting
fig, (ax2, ax3, ax1) = plt.subplots(1, 3, sharey=True, figsize=(5,4))
ax2.spines['right'].set_visible(False)
ax1.tick_params(axis='y', which='both', bottom=False)
ax1.spines['left'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(axis='y', which='both', left=False, right=False)       
ax2.tick_params(axis='y', which='both', left=False, right=False)   
ax3.tick_params(axis='y', which='both', left=False, right=False)   
ax2.tick_params(axis='x', which='both', top=False)    
ax1.tick_params(axis='x', which='both', top=False)
ax3.tick_params(axis='x', which='both', top=False)


ax2.set_xlim(0, 70)
ax2.set_xticks(np.arange(0, 71, 20))
ax3.set_xlim(1000, 4000)
ax3.set_xticks(np.arange(1000, 4001, 1000))
ax1.set_xlim(20000, 91000)
ax1.set_xticks(np.arange(20000, 91001, 20000))
# bars1 = XX.plot(ax=ax1, kind='barh')
# bars2 = XX.plot(ax=ax2, kind='barh')
sns.barplot(ax=ax2,
    data=XX, x="Sample", y="Marker", hue="Group", ci='se', 
    palette=['black', '#fa164f'], capsize=.2)

sns.barplot(ax=ax1,
    data=XX, x="Sample", y="Marker", hue="Group", ci='se', 
    palette=['black', '#fa164f'], capsize=.2)

sns.barplot(ax=ax3,
    data=XX, x="Sample", y="Marker", hue="Group", ci='se', 
    palette=['black', '#fa164f'], capsize=.2)

ax2.legend_.remove()
ax1.legend_.remove()
ax3.legend_.remove()
ax1.set_ylabel('') 
ax2.set_ylabel('')
ax1.set_xlabel('') 
ax3.set_xlabel('') 
ax2.axes.get_yaxis().set_visible(False)
ax3.axes.get_yaxis().set_visible(False)

for tick in ax1.get_xticklabels():
    tick.set_rotation(-90)

for tick in ax2.get_xticklabels():
    tick.set_rotation(-90)

for tick in ax3.get_xticklabels():
    tick.set_rotation(-90)

d = .03
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)      
# ax1.plot((-d, +d),(1 - d, 1 + d), **kwargs)

kwargs.update(transform=ax3.transAxes, color='k', clip_on=False)  
ax3.plot((1 - d, 1 + d), (-d, +d), **kwargs)  
ax3.plot((-d, +d), (-d, +d), **kwargs)   
# ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

kwargs.update(transform=ax2.transAxes)  
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  
# ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

# plt.tight_layout()
plt.subplots_adjust(left=0.125+0.1, bottom=0.1+0.12, right=0.9+0.05)
plt.rc('font', size=15)
plt.rc('legend', fontsize=15) 
plt.rc('xtick', labelsize=15)
ax2.set_xlabel('% Change PHF-tau$\mathregular{^{high}}$\nvs. PHF-tau$\mathregular{^{low}}$ in Postsynaptic', fontsize=13) 
plt.savefig('figures/barplots/hu_post.pdf')




###################
# For Mouse data
###################

# import data ------------------------------------------------------------
df = pd.read_csv('../raw_data/Toms_mean_intensity_change2.csv')

# change data type a bit
for i in range(1, 10):
    df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)] = df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)].apply(lambda x: re.sub(',', '', x))
    df.loc[:, 'Sample'+str(i)] = df.loc[:, 'Sample'+str(i)].astype('float')

# do human pre
df_human_pre = df.loc[(df.HuMu=='Mu') & (df.Region=='Hippo'), :]
df_human_pre = pd.wide_to_long(df_human_pre, ['Sample'], i=['Marker', 'PrePost', 'Region', 'Group', 'HuMu', 'Compare'], \
    j='sample').reset_index()

# Define the sorter
sorter = ['Ab40', 'Ab42', 'CD56', 'SNAP25', 'CD47', 'VGLUT', 'GAD65', 'GATM']
# Create the dictionary that defines the order for sorting
sorterIndex = dict(zip(sorter, range(len(sorter))))
df_human_pre['Marker_rank'] = df_human_pre['Marker'].map(sorterIndex)
XX = df_human_pre.dropna().sort_values(['Marker_rank', 'Compare', 'Marker'], ascending=[True, True, True])

# start plotting
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(5,4))
ax2.yaxis.set_tick_params(labelright=True)
ax1.yaxis.set_tick_params(labelleft=False)
ax1.spines['right'].set_visible(False)
ax1.tick_params(axis='y', which='both', bottom=False)
ax2.spines['left'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(axis='y', which='both', left=False, right=False)       
ax2.tick_params(axis='y', which='both', left=False, right=False)   
ax2.tick_params(axis='x', which='both', top=False) 
ax1.tick_params(axis='x', which='both', top=False)    


ax2.set_xlim(0,60)
ax2.set_xticks(np.arange(0, 61, 20))
ax1.set_xlim(2500, 6000)
ax1.set_xticks(np.arange(2500, 6000, 1000))
# bars1 = XX.plot(ax=ax1, kind='barh')
# bars2 = XX.plot(ax=ax2, kind='barh')


sns.barplot(ax=ax2,
    data=XX, x="Sample", y="Marker", hue="Compare", ci='se', 
    palette=['#738290', '#C2D8B9'], capsize=.2)

sns.barplot(ax=ax1,
    data=XX, x="Sample", y="Marker", hue="Compare", ci='se', 
    palette=['#738290', '#C2D8B9'], capsize=.2)


ax2.legend_.remove()
ax1.legend().set_title('')
ax1.set_ylabel('')
ax2.set_ylabel('')
ax1.set_xlabel('') 

for tick in ax1.get_xticklabels():
    tick.set_rotation(90)

for tick in ax2.get_xticklabels():
    tick.set_rotation(90)

for tick in ax2.get_yticklabels():
    tick.set_ha('center')

# apply offset transform to all x ticklabels.
dx = 0.45; dy = 0. 
offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
for label in ax2.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)


d = .03
kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
ax2.plot((-d, +d), (-d, +d), **kwargs)      
# ax1.plot((-d, +d),(1 - d, 1 + d), **kwargs)

kwargs.update(transform=ax1.transAxes)  
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  
# ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

ax1.invert_xaxis()
ax2.invert_xaxis()
ax1.legend(loc='lower left')

# plt.tight_layout()
plt.subplots_adjust(right=0.9-0.1, bottom=0.1+0.12, left=0.125-0.05)
plt.rc('font', size=15)
plt.rc('legend', fontsize=15) 
plt.rc('xtick', labelsize=15)
# ax2.set_xlabel('% Change Ab40$\mathregular{^{high}}$/Ab42$\mathregular{^{high}}$\nvs. Ab40$\mathregular{^{low}}$/Ab42$\mathregular{^{low}}$ in Hipp', fontsize=13) 
ax2.set_xlabel('% Change Ab$\mathregular{^{high}}$ vs.\nAb$\mathregular{^{low}}$ in Hipp', fontsize=13)
plt.savefig('figures/barplots/mu_hippo.pdf')



# for right hand side
# import data ------------------------------------------------------------
df = pd.read_csv('../raw_data/Toms_mean_intensity_change2.csv')

# change data type a bit
for i in range(1, 10):
    df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)] = df.loc[~df.loc[:, 'Sample'+str(i)].isna(), 'Sample'+str(i)].apply(lambda x: re.sub(',', '', x))
    df.loc[:, 'Sample'+str(i)] = df.loc[:, 'Sample'+str(i)].astype('float')

# do human pre
df_human_pre = df.loc[(df.HuMu=='Mu') & (df.Region=='CC'), :]
# df_human_pre = df_human_pre.dropna(axis=1)
df_human_pre = pd.wide_to_long(df_human_pre, ['Sample'], i=['Marker', 'PrePost', 'Region', 'Group', 'HuMu', 'Compare'], \
    j='sample').reset_index()

# Define the sorter
sorter = ['Ab40', 'Ab42', 'CD56', 'SNAP25', 'CD47', 'VGLUT', 'GAD65', 'GATM']
# Create the dictionary that defines the order for sorting
sorterIndex = dict(zip(sorter, range(len(sorter))))
df_human_pre['Marker_rank'] = df_human_pre['Marker'].map(sorterIndex)
XX = df_human_pre.dropna().sort_values(['Marker_rank', 'Compare', 'Marker'], ascending=[True, True, True])

# start plotting
fig, (ax2, ax1) = plt.subplots(1, 2, sharey=True, figsize=(5,4))
ax2.spines['right'].set_visible(False)
ax1.tick_params(axis='y', which='both', bottom=False)
ax1.spines['left'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(axis='y', which='both', left=False, right=False)       
ax2.tick_params(axis='y', which='both', left=False, right=False)   
ax2.tick_params(axis='x', which='both', top=False)    
ax1.tick_params(axis='x', which='both', top=False)    


ax2.set_xlim(0, 70)
ax2.set_xticks(np.arange(0, 71, 20))
ax1.set_xlim(1000, 6000)
ax1.set_xticks(np.arange(1000, 6000, 1000))
# bars1 = XX.plot(ax=ax1, kind='barh')
# bars2 = XX.plot(ax=ax2, kind='barh')
sns.barplot(ax=ax2,
    data=XX, x="Sample", y="Marker", hue="Compare", ci='se', 
    palette=['#738290', '#C2D8B9'], capsize=.2)


sns.barplot(ax=ax1,
    data=XX, x="Sample", y="Marker", hue="Compare", ci='se', 
    palette=['#738290', '#C2D8B9'], capsize=.2)

ax2.legend_.remove()
ax1.legend_.remove()
ax1.set_ylabel('') 
ax2.set_ylabel('')
ax1.set_xlabel('') 
ax2.axes.get_yaxis().set_visible(False)

for tick in ax1.get_xticklabels():
    tick.set_rotation(-90)

for tick in ax2.get_xticklabels():
    tick.set_rotation(-90)

d = .03
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)      
# ax1.plot((-d, +d),(1 - d, 1 + d), **kwargs)


kwargs.update(transform=ax2.transAxes)  
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  
# ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

# plt.tight_layout()
plt.subplots_adjust(left=0.125+0.1, bottom=0.1+0.12, right=0.9+0.05)
plt.rc('font', size=15)
plt.rc('legend', fontsize=15)
plt.rc('xtick', labelsize=15)
ax2.set_xlabel('% Change Ab$\mathregular{^{high}}$ vs.\nAb$\mathregular{^{low}}$ in Cerebral Cortex', fontsize=13)
plt.savefig('figures/barplots/mu_CC.pdf')
