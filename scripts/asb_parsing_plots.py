import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact
import pyranges as pr


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

# combine peptide and minigene?
# all_detail = pd.concat([pept_detail, mini_detail], axis=0)

# split out TCR beta info
pept_detail[['CDR3-BETA', 'TCRBV', 'TCRBJ']] = pept_detail['TCR BioIdentity'].str.split('+', expand=True)
# unique CDR3 seq
pept_uniq_idx = ~pept_detail['CDR3-BETA'].duplicated()

# pept_detail = pd.read_csv('../../Desktop/COVID19/Data/data_parsed.csv')

# load mcpas data
mcpas = pd.read_csv('McPAS-TCR.csv', low_memory=False, encoding='iso-8859-1')
# filter for human with missing data
mcpas = mcpas.loc[(mcpas['Species'] == 'Human') & ~mcpas['CDR3.beta.aa'].isna() & ~mcpas['Epitope.peptide'].isna() & (mcpas['Category'] == 'Pathogens'), ]
# remove duplicates
mcpas = mcpas.loc[~mcpas[['CDR3.beta.aa', 'Epitope.peptide']].duplicated(), :]

# count unique mcpas records that were found in ImmunoCode COVID screen
mcpas_counts = pd.concat([mcpas.groupby(['Pathology', 'Antigen.protein', 'Epitope.peptide']).size(),
                          mcpas.loc[mcpas['CDR3.beta.aa'].isin(pept_detail['CDR3-BETA']), :].groupby(['Pathology', 'Antigen.protein', 'Epitope.peptide']).size()],
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
_, ax = plt.subplots(ncols=3)
yticks = np.arange(mcpas_counts.shape[0])
yticklabels = mcpas_counts['Pathology'] + ' (' + mcpas_counts['Antigen.protein'] + '): ' + mcpas_counts['Epitope.peptide']
ax[0].barh(yticks, mcpas_counts['mcpas_baseline_prop'], height=0.5)
ax[0].set(yticks=yticks, yticklabels=yticklabels, xlim=[0, 0.4])
ax[1].barh(yticks, mcpas_counts['mcpas_covid_prop'], height=0.5)
ax[1].set(yticks=yticks, yticklabels='', xlim=[0, 0.4])
ax[2].barh(yticks, mcpas_counts['mcpas_delta'], height=0.5)
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
pept_covid_ord_uniq = ~pept_detail[['CDR3-BETA', 'ORF']].duplicated()

covid_orf_counts = pd.concat([pept_detail.loc[pept_covid_ord_uniq, :].groupby('ORF').size(),
                              pept_detail.loc[pept_covid_ord_uniq & pept_detail['CDR3-BETA'].isin(mcpas.loc[mcpas['Epitope.peptide'] == pept_idx, 'CDR3.beta.aa']), :].groupby('ORF').size()],
                             axis=1).fillna(0).reset_index()
covid_orf_counts.columns = ['ORF', 'orfs_baseline', 'orfs_%s' % pept_idx]

covid_orf_counts['orfs_baseline_prop'] = covid_orf_counts['orfs_baseline'] / covid_orf_counts['orfs_baseline'].sum()
covid_orf_counts['orfs_%s_prop' % pept_idx] = covid_orf_counts['orfs_%s' % pept_idx] / covid_orf_counts['orfs_%s' % pept_idx].sum()
covid_orf_counts['delta'] = covid_orf_counts['orfs_%s_prop' % pept_idx] - covid_orf_counts['orfs_baseline_prop']
covid_orf_counts = covid_orf_counts.sort_values('delta')

# plot
_, ax = plt.subplots(ncols=3)
yticks = np.arange(covid_orf_counts.shape[0])
yticklabels = covid_orf_counts.index.values
ax[0].barh(yticks, covid_orf_counts['orfs_baseline_prop'], height=0.5)
ax[0].set(yticks=yticks, yticklabels=yticklabels, xlim=[0, 0.4])
ax[1].barh(yticks, covid_orf_counts['orfs_%s_prop' % pept_idx], height=0.5)
ax[1].set(yticks=yticks, yticklabels='', xlim=[0, 0.4])
ax[2].barh(yticks, covid_orf_counts['delta'], height=0.5)
plt.tight_layout()
covid_orf_counts['ORF']

orf_idx = 'surface glycoprotein'
tab_fischer = pd.concat([covid_orf_counts.loc[covid_orf_counts['ORF'] == orf_idx, ['orfs_baseline', 'orfs_%s' % pept_idx]].T,
                         covid_orf_counts.loc[~(covid_orf_counts['ORF'] == orf_idx), ['orfs_baseline', 'orfs_%s' % pept_idx]].sum(axis=0)], axis=1)
tab_fischer.columns = [orf_idx, 'not_%s' % orf_idx]
tab_fischer.loc['orfs_not_%s' % pept_idx] = tab_fischer.loc['orfs_baseline', :] - tab_fischer.loc['orfs_%s' % pept_idx, :]

fisher_exact(tab_fischer.loc[['orfs_%s' % pept_idx, 'orfs_not_%s' % pept_idx], :])
tab_fischer.loc[['orfs_%s' % pept_idx, 'orfs_not_%s' % pept_idx], :] / tab_fischer.loc['orfs_baseline', :]
