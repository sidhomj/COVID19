"""
Contains ancillary methods.
"""

import numpy as np
import colorsys
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
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

def BarPlot(df_agg):
    df_agg.sort_values(by=['orf_name','counts'],inplace=True,ascending = False)
    df_agg.rename(columns = {'orf_name':'ORF'},inplace=True)
    sns.barplot(data=df_agg,x='peptide',y='counts',order=df_agg['peptide'],hue='ORF',dodge=False)
    plt.xticks(rotation=90)
    plt.subplots_adjust(bottom=0.3)
    plt.xlabel('')
    plt.ylabel('')
