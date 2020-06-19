import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact


pept_detail = pd.read_csv('../Data/data_parsed.csv')

# load mcpas data
mcpas = pd.read_csv('../Data/McPAS-TCR.csv', low_memory=False, encoding='iso-8859-1')
# filter for human with missing data
mcpas = mcpas.loc[(mcpas['Species'] == 'Human') & ~mcpas['CDR3.beta.aa'].isna() & ~mcpas['Epitope.peptide'].isna() & (mcpas['Category'] == 'Pathogens'), ]
# remove duplicates
mcpas = mcpas.loc[~mcpas[['CDR3.beta.aa', 'Epitope.peptide']].duplicated(), :]

total_intersection = len(np.intersect1d(mcpas['CDR3.beta.aa'],pept_detail['beta_sequences']))

# count unique mcpas records that were found in ImmunoCode COVID screen
mcpas_counts = pd.concat([mcpas.groupby(['Pathology', 'Antigen.protein', 'Epitope.peptide']).size(),
                          mcpas.loc[mcpas['CDR3.beta.aa'].isin(pept_detail['beta_sequences']), :].groupby(['Pathology', 'Antigen.protein', 'Epitope.peptide']).size()],
                         axis=1).fillna(0).rename(columns={0: 'mcpas_baseline', 1:'mcpas_covid'})

# limit to records seen in ImmuneCode COVID screen
mcpas_counts = mcpas_counts.loc[mcpas_counts['mcpas_covid'] > 0, :]

# get proportions
mcpas_counts['mcpas_baseline_prop'] = mcpas_counts['mcpas_baseline'] / mcpas_counts['mcpas_baseline'].sum()
mcpas_counts['mcpas_covid_prop'] = mcpas_counts['mcpas_covid'] / mcpas_counts['mcpas_covid'].sum()
mcpas_counts['mcpas_delta'] = mcpas_counts['mcpas_covid_prop'] - mcpas_counts['mcpas_baseline_prop']

# sort and release multiindex
mcpas_counts = mcpas_counts.sort_values('mcpas_delta').reset_index()

# plot
_, ax = plt.subplots(ncols=3,figsize=(12,10))
yticks = np.arange(mcpas_counts.shape[0])
yticklabels = mcpas_counts['Pathology'] + ' (' + mcpas_counts['Antigen.protein'] + ')\n ' + mcpas_counts['Epitope.peptide']
ax[0].barh(yticks, mcpas_counts['mcpas_baseline'], height=0.5,color='white')
temp = ax[0].twiny()
temp.barh(yticks, mcpas_counts['mcpas_baseline_prop'], height=0.5,color='grey')
temp.set(yticks=yticks, xlim=[0, 0.4])
temp.set_yticklabels(yticklabels,wrap=True)
ax[1].barh(yticks, mcpas_counts['mcpas_covid'], height=0.5,color='white')
temp = ax[1].twiny()
temp.barh(yticks, mcpas_counts['mcpas_covid_prop'], height=0.5,color='grey')
temp.set(yticks=yticks, yticklabels='', xlim=[0, 0.4])
ax[2].barh(yticks, mcpas_counts['mcpas_delta'], height=0.5,color='grey')
ax[2].axvline(0,color='k')
ax[2].set(yticks=yticks, yticklabels='', xlim=[-0.16, 0.16])
ax[2].xaxis.tick_top()
plt.tight_layout()

pept_idx = 'GILGFVFTL'
tab_fischer = pd.concat([mcpas_counts.loc[mcpas_counts['Epitope.peptide'] == pept_idx, ['mcpas_baseline', 'mcpas_covid']].T,
                         mcpas_counts.loc[~(mcpas_counts['Epitope.peptide'] == pept_idx), ['mcpas_baseline', 'mcpas_covid']].sum(axis=0)], axis=1)
tab_fischer.columns = [pept_idx, 'not_%s' % pept_idx]
tab_fischer.loc['mcpas_not_covid'] = tab_fischer.loc['mcpas_baseline', :] - tab_fischer.loc['mcpas_covid', :]
fisher_exact(tab_fischer.loc[['mcpas_covid', 'mcpas_not_covid'], :])
tab_fischer.loc[['mcpas_covid', 'mcpas_not_covid'], :] / tab_fischer.loc['mcpas_baseline', :]

## distribution of ORFs for GILGFVFTL

# unique CDR3 ORF records
pept_covid_orf_uniq = ~pept_detail[['beta_sequences', 'orf_name']].duplicated()

covid_orf_counts = pd.concat([pept_detail.loc[pept_covid_orf_uniq, :].groupby('orf_name').size(),
                              pept_detail.loc[pept_covid_orf_uniq & pept_detail['beta_sequences'].isin(mcpas['CDR3.beta.aa']), :].groupby('orf_name').size()],
                             axis=1).fillna(0).reset_index()

covid_orf_counts.columns = ['orf_name', 'orfs_baseline', 'orfs']

covid_orf_counts['orfs_baseline_prop'] = covid_orf_counts['orfs_baseline'] / covid_orf_counts['orfs_baseline'].sum()
covid_orf_counts['orfs_prop'] = covid_orf_counts['orfs'] / covid_orf_counts['orfs'].sum()
covid_orf_counts['delta'] = covid_orf_counts['orfs_prop'] - covid_orf_counts['orfs_baseline_prop']
covid_orf_counts = covid_orf_counts.sort_values('delta')

_, ax = plt.subplots(ncols=3,figsize=(12,10))
yticks = np.arange(covid_orf_counts.shape[0])
yticklabels = covid_orf_counts['orf_name']
ax[0].barh(yticks, covid_orf_counts['orfs_baseline'], height=0.5,color='white')
temp = ax[0].twiny()
temp.barh(yticks, covid_orf_counts['orfs_baseline_prop'], height=0.5,color='grey')
temp.set(yticks=yticks, xlim=[0, 0.5])
temp.set_yticklabels(yticklabels,wrap=True)
ax[1].barh(yticks, covid_orf_counts['orfs'], height=0.5,color='white')
temp = ax[1].twiny()
temp.barh(yticks, covid_orf_counts['orfs_prop'], height=0.5,color='grey')
temp.set(yticks=yticks, yticklabels='', xlim=[0, 0.5])
ax[2].barh(yticks, covid_orf_counts['delta'], height=0.5,color='grey')
ax[2].axvline(0,color='k')
ax[2].set(yticks=yticks, yticklabels='', xlim=[-0.2, 0.2])
ax[2].xaxis.tick_top()
plt.tight_layout()

orf_idx = 'surface glycoprotein'
tab_fischer = pd.concat([covid_orf_counts.loc[covid_orf_counts['orf_name'] == orf_idx, ['orfs_baseline', 'orfs']].T,
                         covid_orf_counts.loc[~(covid_orf_counts['orf_name'] == orf_idx), ['orfs_baseline', 'orfs']].sum(axis=0)], axis=1)
tab_fischer.columns = [orf_idx, 'not']
tab_fischer.loc['orfs_not'] = tab_fischer.loc['orfs_baseline', :] - tab_fischer.loc['orfs', :]
fisher_exact(tab_fischer.loc[['orfs', 'orfs_not'], :])
tab_fischer.loc[['orfs', 'orfs_not'], :] / tab_fischer.loc['orfs_baseline', :]
