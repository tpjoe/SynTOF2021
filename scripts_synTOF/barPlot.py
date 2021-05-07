# import lib
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('classic')
import numpy as np
import re
import seaborn as sns
# pd.options.display.max_rows = 999
import matplotlib

# import data ------------------------------------------------------------
df = pd.read_csv('../raw_data/Toms_HuMu_Compare.csv')

# do human
df_human_pre = df.loc[(df.HuMu=='Hu'), :]
df_human_pre = pd.wide_to_long(df_human_pre, ['Mean'], i=['Marker', 'Region', 'HuMu', ], j='sample').reset_index()

df_human_pre.loc[~df_human_pre.Marker.isin(['PHF-tau', 'DJ1']), 'Mean'] = np.nan

# rename vglut and region
df_human_pre.loc[df_human_pre.Marker=='VGLUT', 'Marker'] = 'vGLUT'
df_human_pre.Region = df_human_pre.Region.map({'BA9':'BA9', 'Hipp':'Hippocampus'})

# remove AS and LRRK2
df_human_pre = df_human_pre.loc[df_human_pre.Marker!='LRRK2', :]
df_human_pre = df_human_pre.loc[df_human_pre.Marker!='AS', :]

axis_color = ['#727273']*13 + ['#d95743']*5 + ['#3251a1']*5 + ['#579660']*6 + ['#6f4a94']*3



# start plotting
fig, ax1 = plt.subplots(1, 1, sharey=True, figsize=(4,11))
ax1.yaxis.set_tick_params(labelleft=False, labelright=True)
ax1.spines['left'].set_visible(False)
ax1.tick_params(axis='y', which='both', bottom=False)
# ax1.spines['left'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.tick_params(axis='y', which='both', left=False, right=True)       
ax1.tick_params(axis='x', which='both', top=False)    

sns.barplot(ax=ax1,
    data=df_human_pre, x="Mean", y="Marker", hue="Region", ci='se', 
    palette=['#513B56', '#348AA7'], capsize=.1)

ax1.legend(loc='top left').set_title('')
ax1.set_xlabel('% Change ADNC vs. Control\n(Average $\pm$ SEM)', fontsize=13) 
ax1.set_ylabel('') 

for tick in ax1.get_xticklabels():
    tick.set_rotation(90)

for i, tick in enumerate(ax1.get_yticklabels()):
    tick.set_ha('center')
    tick.set_color(axis_color[i])

# apply offset transform to all x ticklabels.
dx = 0.7; dy = 0. 
offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
for label in ax1.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

ax1.invert_xaxis()
plt.vlines(x=0, ymin=0, ymax=32, linestyles='--')
# ax2.invert_xaxis()

# plt.tight_layout()
plt.subplots_adjust(right=0.9-0.26, bottom=0.1+0.02, top=0.9+0.07, left=0.125-0.07)
plt.rc('font', size=13)
plt.rc('legend', fontsize=13) 
plt.rc('xtick', labelsize=13)
plt.savefig('figures/barplots/Human.pdf')



#---------------------------------------------------------------------------

# do mouse
df_mouse = df.loc[(df.HuMu=='Mu'), :]
df_mouse = pd.wide_to_long(df_mouse, ['Mean'], i=['Marker', 'Region', 'HuMu', ], j='sample').reset_index()


# rename vglut and region
df_mouse.loc[df_mouse.Marker=='VGLUT', 'Marker'] = 'vGLUT'
df_mouse.Region = df_mouse.Region.map({'CC':'CC', 'Hipp':'Hippocampus'})

# remove AS and LRRK2
df_mouse = df_mouse.loc[df_mouse.Marker!='LRRK2', :]
df_mouse = df_mouse.loc[df_mouse.Marker!='AS', :]
df_mouse = df_mouse.sort_values(['Region'])

axis_color = ['#727273']*13 + ['#d95743']*5 + ['#3251a1']*5 + ['#579660']*6 + ['#6f4a94']*3

# start plotting
fig, ax1 = plt.subplots(1, 1, sharey=True, figsize=(4,11))
ax1.yaxis.set_tick_params(labelleft=False, labelright=True)
ax1.spines['right'].set_visible(False)
ax1.tick_params(axis='y', which='both', bottom=False)
# ax1.spines['left'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.tick_params(axis='y', which='both', left=True, right=False)       
ax1.tick_params(axis='x', which='both', top=False)    

sns.barplot(ax=ax1,
    data=df_mouse, x="Mean", y="Marker", hue="Region", ci='se', 
    palette=['#513B56', '#348AA7'], capsize=.1)

ax1.set_xlabel('% Change PS/APP vs. WT\n(Average $\pm$ SEM)', fontsize=13) 


for tick in ax1.get_xticklabels():
    tick.set_rotation(-90)

# apply offset transform to all x ticklabels.
dx = 0.7; dy = 0. 
offset = matplotlib.transforms.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
for label in ax1.yaxis.get_majorticklabels():
    label.set_transform(label.get_transform() + offset)

# ax1.invert_xaxis()
plt.vlines(x=0, ymin=0, ymax=32, linestyles='--')
ax1.set(xscale="log")
# ax2.invert_xaxis()
ax1.legend(loc='center right').set_title('')

L=ax1.legend()
L.get_texts()[0].set_text('Cerebral\nCortex')
ax1.set_ylabel('') 


# plt.tight_layout()
plt.subplots_adjust(right=0.9+0.07, bottom=0.1+0.02, top=0.9+0.07, left=0.125+0.26)
plt.rc('font', size=13)
plt.rc('legend', fontsize=13) 
plt.rc('xtick', labelsize=13)
plt.savefig('figures/barplots/Mouse.pdf')




