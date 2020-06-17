import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


# data files
subj = pd.read_csv('subject-metadata.csv')
pept_hits = pd.read_csv('peptide-hits.csv')
pept_detail = pd.read_csv('peptide-detail.csv')
mini_hits = pd.read_csv('minigene-hits.csv')
mini_detail = pd.read_csv('minigene-detail.csv')

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
mcpas = pd.read_csv('McPAS-TCR.csv', low_memory=False, encoding='iso-8859-1')
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
_, ax = plt.subplots(ncols=2)
yticks = np.arange(tab.shape[0])
yticklabels = tab['Pathology'] + ' (' + tab['Antigen.protein'] + ')'
ax[0].barh(yticks, tab['all mcpas'].values, height=0.5)
ax[0].set(yticks=yticks, yticklabels=yticklabels, title='All TCR in McPAS')
ax[1].barh(np.arange(tab.shape[0]), tab['covid specific'].values, height=0.5)
ax[1].set(yticks=yticks, yticklabels='', title='COVID Specific TCR in McPAS')
plt.tight_layout()

pd.concat([tab.loc[(tab['Pathology'] == 'Influenza') & (tab['Antigen.protein'] == 'Matrix protein (M1)'), ['covid specific', 'all mcpas']].T,
           tab.loc[~((tab['Pathology'] == 'Influenza') & (tab['Antigen.protein'] == 'Matrix protein (M1)')), ['covid specific', 'all mcpas']].sum(axis=0)], axis=1)

