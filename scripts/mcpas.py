import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from utils import *
import matplotlib
matplotlib.rc('font', family='Arial')

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
label_dict = dict(zip(data['beta_sequences'],data['orf_name']))
labels = list(map(label_dict.get,data_piv.index))

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

DFs = []
for cl in np.unique(idx):
    tcr = data_piv.index[idx == cl]
    DFs.append(data[data['beta_sequences'].isin(tcr)])

cluster_sel = 1
df_agg = DFs[cluster_sel].groupby(['peptide', 'orf_name']).agg({'counts': 'sum'}).reset_index()
BarPlot(df_agg)

cluster_sel = 2
df_agg = DFs[cluster_sel].groupby(['peptide', 'orf_name']).agg({'counts': 'sum'}).reset_index()
BarPlot(df_agg)









