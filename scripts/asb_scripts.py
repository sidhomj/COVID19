import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact
from scripts.utils import delta_bar_plots


pept_detail = pd.read_csv('./Data/data_parsed.csv')

# load mcpas data
mcpas = pd.read_csv('./Data/McPAS-TCR.csv', low_memory=False, encoding='iso-8859-1')
# filter for human with missing data
mcpas = mcpas.loc[(mcpas['Species'] == 'Human') & ~mcpas['CDR3.beta.aa'].isna() & ~mcpas['Epitope.peptide'].isna() & (mcpas['Category'] == 'Pathogens'), ]
# remove duplicates
mcpas = mcpas.loc[~mcpas[['CDR3.beta.aa', 'Epitope.peptide']].duplicated(), :]

# number of unique tcr common to both mcpas and covid
total_intersection = len(np.intersect1d(mcpas['CDR3.beta.aa'],pept_detail['beta_sequences']))

# count unique mcpas records that were found in ImmunoCode COVID screen
mcpas_covid = pd.concat([mcpas.groupby(['Pathology', 'Antigen.protein', 'Epitope.peptide']).size(),
                         mcpas.loc[mcpas['CDR3.beta.aa'].isin(pept_detail['beta_sequences']), :].groupby(['Pathology', 'Antigen.protein', 'Epitope.peptide']).size()],
                        axis=1).fillna(0).rename(columns={0: 'baseline', 1: 'covid'})

# limit to records seen in ImmuneCode COVID screen
mcpas_covid = mcpas_covid.loc[mcpas_covid['covid'] > 0, :]
mcpas_covid['not_covid'] = mcpas_covid['baseline'] - mcpas_covid['covid']

# get proportions
mcpas_covid['baseline_prop'] = mcpas_covid['baseline'] / mcpas_covid['baseline'].sum()
mcpas_covid['covid_prop'] = mcpas_covid['covid'] / mcpas_covid['covid'].sum()
mcpas_covid['delta'] = mcpas_covid['covid_prop'] - mcpas_covid['baseline_prop']

# sort and release multiindex
mcpas_covid = mcpas_covid.sort_values('delta').reset_index()

# fishers
mcpas_covid['fisher_p'] = None
for idx in mcpas_covid.index:
    idx_bool = mcpas_covid.index == idx
    mcpas_covid.loc[idx, 'fisher_p'] = fisher_exact(np.stack([mcpas_covid.loc[idx_bool, ['covid', 'not_covid']].sum(axis=0).values,
                                                              mcpas_covid.loc[~idx_bool, ['covid', 'not_covid']].sum(axis=0).values], axis=1))[1]

# plot
delta_bar_plots(baseline=mcpas_covid[['baseline', 'baseline_prop']].values,
                signal=mcpas_covid[['covid', 'covid_prop']].values,
                yticklabels = mcpas_covid['Pathology'] + ' (' + mcpas_covid['Antigen.protein'] + ')\n ' + mcpas_covid['Epitope.peptide'],
                max_proporption=.4, max_delta=0.18)
# resize then call
plt.tight_layout()

# unique CDR3 ORF records
pept_covid_orf_uniq = ~pept_detail[['beta_sequences', 'orf_name']].duplicated()

pept_idx = 'GILGFVFTL'
covid_orf_counts = pd.concat([pept_detail.loc[pept_covid_orf_uniq, :].groupby('orf_name').size(),
                              pept_detail.loc[pept_covid_orf_uniq & pept_detail['beta_sequences'].isin(mcpas.loc[mcpas['Epitope.peptide'] == pept_idx, 'CDR3.beta.aa']), :].groupby('orf_name').size()],
                             axis=1).fillna(0).reset_index()

covid_orf_counts.columns = ['orf_name', 'baseline', 'orf']
covid_orf_counts['not_orf'] = covid_orf_counts['baseline'] - covid_orf_counts['orf']
covid_orf_counts['baseline_prop'] = covid_orf_counts['baseline'] / covid_orf_counts['baseline'].sum()
covid_orf_counts['orf_prop'] = covid_orf_counts['orf'] / covid_orf_counts['orf'].sum()
covid_orf_counts['delta'] = covid_orf_counts['orf_prop'] - covid_orf_counts['baseline_prop']
covid_orf_counts = covid_orf_counts.sort_values('delta')

# fishers
covid_orf_counts['fisher_p'] = None
for idx in covid_orf_counts.index:
    idx_bool = covid_orf_counts.index == idx
    covid_orf_counts.loc[idx, 'fisher_p'] = fisher_exact(np.stack([covid_orf_counts.loc[idx_bool, ['orf', 'not_orf']].sum(axis=0).values,
                                                                   covid_orf_counts.loc[~idx_bool, ['orf', 'not_orf']].sum(axis=0).values], axis=1))[1]

# plot
delta_bar_plots(baseline=covid_orf_counts[['baseline', 'baseline_prop']].values,
                signal=covid_orf_counts[['orf', 'orf_prop']].values,
                yticklabels = covid_orf_counts['orf_name'],
                max_proporption=.7, max_delta=0.6)
# resize then call
plt.tight_layout()
