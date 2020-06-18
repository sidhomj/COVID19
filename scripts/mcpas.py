"""
Code to generate results and figures from manuscript.
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

# data files
subj = pd.read_csv('../Data/subject-metadata.csv')
pept_hits = pd.read_csv('../Data/peptide-hits.csv')
pept_detail = pd.read_csv('../Data/peptide-detail.csv')
mini_hits = pd.read_csv('../Data/minigene-hits.csv')
mini_detail = pd.read_csv('../Data/minigene-detail.csv')

# columns name adjustments
pept_detail.rename(columns={'ORF Coverage': 'ORF', 'Amino Acids': 'Amino Acid'}, inplace=True)
mini_detail.drop(columns=['ORF Genebank ID'], inplace=True)
(pept_detail.columns == mini_detail.columns).all()
pept_detail['assay type'] = 'peptide'
mini_detail['assay type'] = 'minigene'

# combine peptide and minigene data and split out TCR beta info
all_detail = pd.concat([pept_detail, mini_detail], axis=0)
all_detail[['CDR3-BETA', 'TCRBV', 'TCRBJ']] = all_detail['TCR BioIdentity'].str.split('+', expand=True)

# load mcpas data and filter
mcpas = pd.read_csv('../Data/McPAS-TCR.csv', low_memory=False, encoding='iso-8859-1')
idx = (mcpas['Species'] == 'Human') & ~mcpas['CDR3.beta.aa'].isna() & (mcpas['Category'] == 'Pathogens')
mcpas_human = mcpas.loc[idx, ]

# background distribution of data in mcpas
mcpas_bg = mcpas_human.groupby(['Pathology', 'Antigen.protein']).size()

# cross mcpass TCR beta against covid hits (exact sequence  match ...)
mcpas_human_X_covid = mcpas_human['CDR3.beta.aa'].values[:, np.newaxis] == all_detail['CDR3-BETA'].values[np.newaxis, :]
# covid specific distribution of data in mcpas
mcpas_covid = mcpas_human.loc[mcpas_human_X_covid.any(axis=1), :].groupby(['Pathology', 'Antigen.protein']).size()

# combine tables
tab = pd.concat([mcpas_bg, mcpas_covid], join='outer', axis=1).rename(columns={0: 'all mcpas', 1: 'covid specific'}).fillna(0).reset_index()
tab = tab.sort_values(['Pathology', 'all mcpas'], ascending=[True, False])

# plot
fig, ax = plt.subplots(ncols=2,figsize=(10,8))
yticks = np.arange(tab.shape[0])
yticklabels = tab['Pathology'] + ' (' + tab['Antigen.protein'] + ')'
ax[0].barh(yticks, tab['all mcpas'].values, height=0.5,color='grey')
ax[0].set(yticks=yticks, yticklabels=yticklabels, title='All TCRs in McPAS')
ax[1].barh(np.arange(tab.shape[0]), tab['covid specific'].values, height=0.5,color='grey')
ax[1].set(yticks=yticks, yticklabels='', title='COVID Specific TCRs in McPAS')
plt.tight_layout()
fig.savefig('../Results/1A.eps')

#pull parsed data
data = pd.read_csv('../Data/data_parsed.csv')
data['counts'] =  1
mcpas = pd.read_csv('../Data/McPAS-TCR.csv')
mcpas = mcpas[mcpas['Species']=='Human']
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











