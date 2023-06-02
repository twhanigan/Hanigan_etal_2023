import numpy as np
import pandas as pd
import math 
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy
from scipy import stats
from scipy.stats import pearsonr
import warnings
from scipy.stats import ttest_ind
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lmfit.models import LinearModel, StepModel
from lmfit import Model
from matplotlib import gridspec
import re
from adjustText import adjust_text
from pylab import *
from glob import glob
import matplotlib.style as mplstyle
import matplotlib.colors as mcolors
import itertools
plt.rcParams['agg.path.chunksize'] = 10000
import matplotlib.font_manager as font_manager

df_subset= pd.read_csv('Merged_Aggregated_Frame_c.csv',)#sep = '\t',comment='@',skip_blank_lines=True)

#df_subset = df_subset.sort_values(['Chromosome','Start']).reset_index(drop=True)
#df_subset["group"]=(df_subset["Start"]>df_subset["End"].shift()).cumsum()
## this returns min value of "START" column from a group and max value fro m "FINISH"
#result=df_subset.groupby(['Chromosome',"group"]).agg({"Start":"min", "End": "max",'Percentage':lambda x : abs(x).max()})
#print(result)
#result.to_csv('Merged_Aggregated_Frame.csv')
#Order the columns
df_subset.loc[(df_subset['Strand']=='-'),'Percentage'] = df_subset.loc[(df_subset['Strand']=='-'),'Percentage']*-1
col_order=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrM']#'chrX','chrY',
sorterIndex = dict(zip(col_order, range(len(col_order))))
df_subset['Tm_Rank'] = df_subset['Chromosome'].map(sorterIndex)
df_subset = df_subset.sort_values(by=['Tm_Rank'],ascending = True).copy()
df_subset.drop(['Tm_Rank'], 1, inplace = True)
df_subset = df_subset.loc[df_subset['Chromosome'].isin(col_order)]
chrome_sizes = pd.read_excel('Chrom_Sizes.xlsx')
#chrome_sizes = chrome_sizes[chrome_sizes['Chromosome']!='chr18']
chrome_sizes = chrome_sizes.loc[chrome_sizes['Chromosome'].isin(col_order)]
widths = chrome_sizes['Size'].to_list()
print(df_subset.Chromosome.unique())
heights = [0.5]
df_subset.loc[(df_subset['Chromosome']=='chr4'),'Percentage'] = df_subset.loc[(df_subset['Chromosome']=='chr4'),'Percentage']/10
#df_subset.loc[(df_subset['Chromosome']=='chrM'),'Percentage'] = df_subset.loc[(df_subset['Chromosome']=='chrM'),'Percentage']*1.5
#df_subset.loc[(df_subset['Chromosome']=='chr16'),'Percentage'] = df_subset.loc[(df_subset['Chromosome']=='chr16'),'Percentage']/2
#df_subset.loc[(df_subset['Chromosome']=='chr1'),'Percentage'] = df_subset.loc[(df_subset['Chromosome']=='chr1'),'Percentage']/1.5
lister= ['chr14','chrM']
df_subset.loc[~(df_subset['Chromosome'].isin(lister)),'Percentage'] = df_subset.loc[~(df_subset['Chromosome'].isin(lister)),'Percentage']/3

#df_subset.loc[(df_subset['Chromosome']=='chrM')].to_csv('merged_frame_recurrent_chrM.csv')
spec = dict(height_ratios=heights,width_ratios=widths, hspace=0.0,wspace=0.05)
fig,axes = plt.subplots(ncols=len(df_subset['Chromosome'].unique()),nrows=1,sharex=False, sharey=False,constrained_layout=False,figsize=(10,3.33),gridspec_kw=spec)#(20,5)(15,5)
colors = ["firebrick",'steelblue']
colors_b = ['white','grey']
color_pal= sns.color_palette(colors, n_colors=len(col_order))
color_pal_b = sns.color_palette(colors, n_colors=2)
color_pal_background = sns.color_palette(colors_b, n_colors=12)
palette = itertools.cycle(color_pal_b)
palette_background = itertools.cycle(colors_b)
n=-1
for group, ax in zip(df_subset['Chromosome'].unique(), axes.flat):
    c_back = next(palette_background)
    maximum = df_subset['End'].loc[(df_subset['Chromosome']==group)].max()
    minimum = df_subset['Start'].loc[(df_subset['Chromosome']==group)].min()
    pos_signal = df_subset.loc[(df_subset['Chromosome']==group)].copy()
    neg_signal = df_subset.loc[(df_subset['Chromosome']==group)].copy()
    pos_signal.loc[(pos_signal['Percentage'] < 0),'Percentage'] = 0
    neg_signal.loc[(neg_signal['Percentage'] > 0),'Percentage'] = 0
    x = np.linspace(maximum,minimum,num=len(pos_signal))
    ax.fill_between(x,pos_signal['Percentage'],0,step='mid',facecolor='red',alpha=0.7)#step='mid',
    ax.plot(x,pos_signal['Percentage'],linewidth=0.02,color='firebrick',)#step='mid',
    ax.fill_between(x,neg_signal['Percentage'],0,step='mid',color='dodgerblue',alpha=0.7)#facecolor=(0.11764705882352941, 0.5647058823529412, 0.95))
    ax.plot(x,neg_signal['Percentage'],linewidth=0.02,color='steelblue',)
    sns.despine(left=True)#sns.stripplot(x='Compound',y='value',hue='dataset',edgecolor='black',linewidth=1,data=concat_subset_c,dodge=True,palette=color_pal,size=6,jitter=0.3,ax=ax,alpha=0.9,hue_order = ['Sensitive(NH4)','Sensitive','Resistant'])
    ax.set(xlabel=None,ylabel=None,facecolor=c_back)
    if ((group == 'chr14')| (group =='chrM')):
        ax.set_xlabel(group,fontsize=24,rotation=40,ha="right",color='dodgerblue',fontweight='bold')
    else:
        ax.set_xlabel(group,fontsize=18,rotation=40,ha="right",color='black',fontweight='bold')
    ax.set_yticks([-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1])
    ax.patch.set_alpha(0.25)
    ax.set_xticks([])
    #ax.set_xlim((minimum,maximum))
    #ax.set_xlabel(group,fontsize=18,rotation=40,ha="right",fontweight='bold')
    ax.axhline(y=0,color='black')
plt.setp(axes, ylim=(-1.1,1.1),yticks=[])
#plt.xticks([])
#axes[1].set_yticks([],fontsize=16,)
axes[0].set_yticks([-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1])
axes[0].set_yticklabels(['1','0.75','0.5','0.25','0','0.25','0.5','0.75','1'],fontsize=16,)

all_axes = fig.get_axes()
#for ax in all_axes[1:]:
#    plt.setp(ax.get_xticklabels(), visible=False)
# show only the outside spines
for ax in all_axes:
    ax.label_outer()
    if ax.is_first_row():
        ax.spines['bottom'].set_visible(True)
#        #ax.spines['bottom'].set_position(('outward',10))
        ax.spines['top'].set_visible(True)
#    if ax.is_last_row():
#        ax.spines['bottom'].set_visible(True)
#        #ax.spines['bottom'].set_position(('outward',10))
    if ax.is_first_col():
        ax.spines['left'].set_visible(True)
#       #ax.spines['left'].set_position(('outward',10))
        #ax.set_yticks([-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1])
        #ax.set_yticklabels(['1','0.75','0.5','0.25','0','0.25','0.5','0.75','1'],fontsize=16,)
    if ax.is_last_col():
        ax.spines['right'].set_visible(True)
    #else:
    #    ax.set_yticks([])

#fig.text(0.5, 0.02, 'Chromosome position (bp)',fontsize=18,fontname="Arial",fontweight='bold', ha='center', va='center',color='black',wrap=True)
fig.text(0.035, 0.5, '% Clones with CNA',fontsize=18,fontname="Arial",fontweight='bold',ha='center', va='center', rotation='vertical',color='black',wrap=True)
plt.savefig('Percentage_Samples_10_op05_Final_H.png',dpi=600,transparent=False,bbox_inches='tight',pad_inches=1)
plt.show()
