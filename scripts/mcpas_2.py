import pandas as pd
import numpy as np
import seaborn as sns
import os
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import colorsys
import umap
from sklearn.cluster import KMeans

def Get_Color_Dict(labels):
    N = len(np.unique(labels))
    HSV_tuples = [(x * 1.0 / N, 1.0, 0.5) for x in range(N)]
    np.random.shuffle(HSV_tuples)
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    color_dict = dict(zip(np.unique(labels), RGB_tuples))
    return color_dict

data = pd.read_csv('Data/data_parsed.csv')
mcpas = pd.read_csv('Data/McPAS-TCR.csv')
mcpas = mcpas[mcpas['Species']=='Human']
path_dict = dict(zip(mcpas['CDR3.beta.aa'],mcpas['Pathology']))
data['Pathology'] = data['beta_sequences'].map(path_dict)
# data = data[data['Cohort'] != 'Healthy (No known exposure)']
# data.dropna(subset=['Pathology'],inplace=True)
pep_dict = dict(zip(mcpas['CDR3.beta.aa'],mcpas['Epitope.peptide']))
data['Path_Peptide'] = data['beta_sequences'].map(pep_dict)
data.dropna(subset=['Path_Peptide'],inplace=True)
data['counts'] =  1

piv_subject = pd.pivot_table(data,index='beta_sequences',columns='Path_Peptide',values='counts',
                             fill_value=0.0,aggfunc='max')
label_dict = dict(zip(data['beta_sequences'],data['orf_name']))
labels = list(map(label_dict.get,piv_subject.index))

# color_dict = Get_Color_Dict(labels)
# row_colors = [color_dict[x] for x in labels]
# CM = sns.clustermap(data=piv_subject.T,cmap='jet')
# ax = CM.ax_heatmap
# ax.set_xticklabels('')
# # ax.set_yticklabels('')
# plt.subplots_adjust(right=0.8)

X = np.array(piv_subject)
idx = KMeans(n_clusters=3).fit_predict(X)
color_dict = Get_Color_Dict(idx)
row_colors = [color_dict[x] for x in idx]

plt.figure()
X_2 = umap.UMAP(random_state=6).fit_transform(np.array(piv_subject))
plt.scatter(X_2[:,0],X_2[:,1],c=row_colors)

CM = sns.clustermap(data=piv_subject.T,cmap='binary',col_colors=row_colors)
ax = CM.ax_heatmap
ax.set_xticklabels('')
plt.subplots_adjust(right=0.8)

DFs = []
for cl in np.unique(idx):
    tcr = piv_subject.index[idx == cl]
    DFs.append(data[data['beta_sequences'].isin(tcr)])

DFs[2]['peptide'].value_counts()
plt.figure()
sns.barplot(data=DFs[2]['orf_name'].value_counts().reset_index(),x='index',y='orf_name',)
plt.xticks(rotation=90)
plt.subplots_adjust(bottom=0.3)

DF = pd.concat(DFs)

plt.figure()
sns.barplot(data=data['orf_name'].value_counts().reset_index(),x='index',y='orf_name')
plt.xticks(rotation=90)
plt.subplots_adjust(bottom=0.3)








