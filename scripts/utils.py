"""
Contains ancillary methods.
"""

import numpy as np
import colorsys
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import logomaker
from matplotlib.ticker import MaxNLocator
matplotlib.rc('font', family='Arial')

def Process_Seq(df,col):
    #Drop null values
    df = df.dropna(subset=[col])

    #strip any white space and remove non-IUPAC characters
    df[col] = df[col].str.strip()
    searchfor = ['\*', 'X', 'O']
    df = df[~df[col].str.contains('|'.join(searchfor))]

    return df

def Get_Color_Dict(labels):
    N = len(np.unique(labels))
    HSV_tuples = [(x * 1.0 / N, 1.0, 0.5) for x in range(N)]
    np.random.shuffle(HSV_tuples)
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    color_dict = dict(zip(np.unique(labels), RGB_tuples))
    return color_dict

def BarPlot(df_agg,figsize=(8,6)):
    plt.figure(figsize=figsize)
    df_agg.sort_values(by=['orf_name','counts'],inplace=True,ascending = False)
    df_agg.rename(columns = {'orf_name':'ORF'},inplace=True)
    ax = sns.barplot(data=df_agg,x='peptide',y='counts',order=df_agg['peptide'],hue='ORF',dodge=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks(rotation=90)
    plt.yticks(fontsize=24)
    plt.subplots_adjust(bottom=0.3)
    plt.xlabel('')
    plt.ylabel('')
    leg = plt.legend(loc='upper right',frameon=False,prop={'size': 16})
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    return leg

def BarPlotCohort(df_agg):
    df_agg.sort_values(by='Cohort', inplace=True, ascending=False)
    plt.figure(figsize=(8,6))
    ax = sns.barplot(data=df_agg, x='Subject', y='counts', hue='Cohort', dodge=False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(rotation=90,fontsize=16)
    plt.yticks(fontsize=24)
    plt.subplots_adjust(bottom=0.2)
    leg = plt.legend(frameon=False,prop={'size': 16},loc='upper right')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    return leg

def Get_Logo_df(motifs_logo,kernel):
    df_motifs = pd.DataFrame(motifs_logo)
    df_motifs = df_motifs[0].apply(lambda x: pd.Series(list(x)))
    df_motifs.fillna(value='X',inplace=True)
    cols = np.unique(df_motifs)
    df_out = pd.DataFrame()
    df_out['pos'] = list(range(kernel))
    for c in cols:
        df_out[c] = None
    df_out.set_index('pos', inplace=True)
    for i in range(kernel):
        temp = df_motifs[i].value_counts()
        for k in np.array(temp.index):
            df_out.loc[i, k] = temp[k] / np.sum(temp)
    df_out.fillna(value=0.0, inplace=True)
    if 'X' in cols:
        df_out.drop(columns=['X'], inplace=True)
    return df_out

def Make_Logo(sel_seq):
    df_logo = Get_Logo_df(sel_seq,np.max([len(x) for x in sel_seq]))
    ax = logomaker.Logo(df_logo,color_scheme='weblogo_protein')
    ax.style_spines(spines=['top', 'right', 'left', 'bottom'], visible=False)
    ax.ax.set_xticks([])
    ax.ax.set_yticks([])
    return ax

def delta_bar_plots(baseline, signal, yticklabels, fig_size=(12, 10), max_proporption=0.4, max_delta=0.18):
    _, ax = plt.subplots(ncols=3, figsize=fig_size)
    yticks = np.arange(len(yticklabels))

    ax[0].barh(yticks, baseline[:, 0], height=0.5, color='white')
    y_twin = ax[0].twiny()
    y_twin.barh(yticks, baseline[:, 1], height=0.5, color='grey')
    y_twin.set(yticks=yticks, yticklabels=yticklabels, xlim=[0, max_proporption])

    ax[1].barh(yticks, signal[:, 0], height=0.5, color='white')
    y_twin = ax[1].twiny()
    y_twin.barh(yticks, signal[:, 1], height=0.5, color='grey')
    y_twin.set(yticks=yticks, yticklabels='', xlim=[0, max_proporption])

    ax[2].barh(yticks, signal[:, 1] - baseline[:, 1], height=0.5, color='grey')
    ax[2].axvline(0, color='k')
    ax[2].set(yticks=yticks, yticklabels='', xlim=[-max_delta, max_delta])
    ax[2].xaxis.tick_top()
