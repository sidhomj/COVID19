"""
Code to generate results and figures from manuscript. This is the same code as is present in the jupyter notebook.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from utils import *
import matplotlib
matplotlib.rc('font', family='Arial')
import logomaker

#pull parsed data
data = pd.read_csv('../Data/data_parsed.csv')
data['counts'] =  1
mcpas = pd.read_csv('../Data/McPAS-TCR.csv')
idx = (mcpas['Species'] == 'Human') & ~mcpas['CDR3.beta.aa'].isna() & (mcpas['Category']=='Pathogens')
mcpas = mcpas.loc[idx, ]
path_dict = dict(zip(mcpas['CDR3.beta.aa'],mcpas['Pathology']))
data['Pathology'] = data['beta_sequences'].map(path_dict)
pep_dict = dict(zip(mcpas['CDR3.beta.aa'],mcpas['Epitope.peptide']))
data['Path_Peptide'] = data['beta_sequences'].map(pep_dict)
data.dropna(subset=['Path_Peptide'],inplace=True)

data_piv = pd.pivot_table(data,index='beta_sequences',columns='Path_Peptide',values='counts',
                             fill_value=0.0,aggfunc='max')

X = np.array(data_piv)
idx = KMeans(n_clusters=3,random_state=2).fit_predict(X)
color_dict = Get_Color_Dict(idx)
col_colors = [color_dict[x] for x in idx]
data_piv_sort = data_piv.iloc[np.flip(np.argsort(idx))]
col_colors = np.array(col_colors)[np.flip(np.argsort(idx))]

CM = sns.clustermap(data=data_piv_sort.T,cmap='binary',col_cluster=False,col_colors=col_colors)
ax = CM.ax_heatmap
CM.cax.set_visible(False)
ax.set_xticklabels('')
plt.subplots_adjust(right=0.8,top=1.1,bottom=0.1)
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticks([])
CM.savefig('../Results/clustermap.eps')

DFs = []
for cl in np.unique(idx):
    tcr = data_piv.index[idx == cl]
    DFs.append(data[data['beta_sequences'].isin(tcr)])

cluster_sel = 1
df_agg = DFs[cluster_sel].groupby(['peptide', 'orf_name']).agg({'counts': 'sum'}).reset_index()
df_agg.sort_values(by=['orf_name','peptide'],inplace=True,ascending=[False,True])
df_agg.reset_index(drop=True,inplace=True)
leg = BarPlot(df_agg)
bbox_to_anchor=[1.15, 1]
leg.set_bbox_to_anchor(bbox_to_anchor)

plt.savefig('../Results/cluster1_pep.eps')
DFs[cluster_sel].sort_values(by=['orf_name','peptide'],inplace=True,ascending=[False,True])
sel = [0,1,2,7]
df_sel = DFs[cluster_sel][DFs[cluster_sel]['peptide'].isin(df_agg['peptide'].iloc[sel])]
sel_seq = np.unique(df_sel['beta_sequences'])
ax = Make_Logo(sel_seq)
ax.fig.savefig('../Results/np177_cluster_1.eps')

df_agg = DFs[cluster_sel].drop_duplicates(subset=['beta_sequences','Subject'])
df_agg = df_agg.groupby(['Subject']).agg({'Cohort':'first','counts':'sum'}).reset_index()
BarPlotCohort(df_agg)
plt.savefig('../Results/cluster1_subj.eps')


cluster_sel = 2
df_agg = DFs[cluster_sel].groupby(['peptide', 'orf_name']).agg({'counts': 'sum'}).reset_index()
df_agg.sort_values(by=['orf_name','peptide'],inplace=True,ascending=[False,True])
BarPlot(df_agg)
plt.savefig('../Results/cluster2_pep.eps')
DFs[cluster_sel].sort_values(by=['orf_name','peptide'],inplace=True,ascending=[False,True])
df_sel = DFs[cluster_sel][DFs[cluster_sel]['peptide'].isin(df_agg['peptide'].iloc[0:5])]
sel_seq = np.unique(df_sel['beta_sequences'])
ax = Make_Logo(sel_seq)
ax.fig.savefig('../Results/m1_cluster_1.eps')

df_agg = DFs[cluster_sel].drop_duplicates(subset=['beta_sequences', 'Subject'])
df_agg = df_agg.groupby(['Subject']).agg({'Cohort': 'first', 'counts': 'sum'}).reset_index()
BarPlotCohort(df_agg)
plt.savefig('../Results/cluster2_subj.eps')











