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

file_names = glob('*.tsv')
seg_names = glob('*.seg')
names = [file.split('.called.seg')[0] for file in seg_names]
sorted_frames =sorted(file_names)#,key=lambda x:x.split('.denoisedCR.tsv')[0]
sorted_segs = sorted(seg_names)#,key=lambda x:x.split('.called.seg')[0]
sorted_names = sorted(names)

#seg_names = [file.split('.denoisedCR.tsv')[0] for file in seg_names]
frames = [pd.read_csv(f, sep = '\t',comment='@',skip_blank_lines=True) for f in sorted_frames]
segs = [pd.read_csv(f, sep = '\t',comment='@',skip_blank_lines=True) for f in sorted_segs]
chromo = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrM']
#for frame in segs:
    #frame.columns = ['CONTIG','START','END','NUM_POINTS_COPY_RATIO','LOG2_COPY_RATIO','CALL',]
    #frame = frame[frame['CONTIG'].isin(chromo)]
    #frame = frame[(frame['CALL']!=0)]
    #frame.loc[(frame['CONTIG']=='chrM'),'NUM_POINTS_COPY_RATIO'] = 1500
segs = [frame.loc[frame['CONTIG'].isin(chromo)] for frame in segs]
segs = [frame[(frame['CALL']!=0)] for frame in segs]
#segs = [frame.loc[(frame['CONTIG']=='+')|(frame['CONTIG']=='-')] for frame in segs]
#segs = [frame.loc[(frame['CONTIG']=='chrM')] for frame in segs]
#segs = [frame.loc[((frame['NUM_POINTS_COPY_RATIO']>100)&((frame['MEAN_LOG2_COPY_RATIO']>0.1)|(frame['MEAN_LOG2_COPY_RATIO']<-0.1)))|(frame['CONTIG']=='chrM')] for frame in segs]
frames = [frame.sort_values(['CONTIG','START']) for frame in frames]#.reset_index(drop=True)
chrome_sizes = pd.read_excel('Chrom_Sizes.xlsx')
chrome_sizes = chrome_sizes.loc[chrome_sizes['Chromosome'].isin(chromo)]
widths = chrome_sizes['Size'].to_list()
#Order the columns
def plot_CNV(df_subset,segments,frame_name):
    col_order=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrM']
    sorterIndex = dict(zip(col_order, range(len(col_order))))
    df_subset['Tm_Rank'] = df_subset['CONTIG'].map(sorterIndex)
    df_subset = df_subset.sort_values(by=['Tm_Rank'],ascending = True)
    df_subset.drop(['Tm_Rank'], 1, inplace = True)
    df_subset.LOG2_COPY_RATIO = df_subset.LOG2_COPY_RATIO.astype(np.float)
    df_subset.START = df_subset.START.astype(np.float)
    heights = [0.5]
    group_num = [0,1]
    df_subset['Group'] = 0
    df_subset.loc[(df_subset['CONTIG']=='chr14'),'Group'] = 1
    df_subset.loc[(df_subset['CONTIG']=='chrM'),'Group'] = 1
    spec = dict(height_ratios=heights,width_ratios=widths, hspace=0.0,wspace=0.05)
    fig,axes = plt.subplots(ncols=len(widths),nrows=1,sharex=False, sharey=False,constrained_layout=False,figsize=(10,3.33),gridspec_kw=spec)#(20,5)(15,5)
    colors = ["grey",'dodgerblue']
    colors_b = ['white','grey']
    color_pal= sns.color_palette(colors, n_colors=len(df_subset['Group'].unique()))
    c1 = sns.color_palette(['dodgerblue','grey'],n_colors=2)
    c2 = sns.color_palette(['grey'],n_colors=1)
    print(df_subset.head())
    color_pal_b = sns.color_palette(colors, n_colors=1)
    color_pal_background = sns.color_palette(colors_b, n_colors=2)
    palette = itertools.cycle(color_pal_b)
    palette_background = itertools.cycle(colors_b)
    n=-1
    for group, ax in zip(df_subset['CONTIG'].unique(), axes.flat):
        c_back = next(palette_background)
        if (group == 'chr14'):
            ax.scatter(x='START',y='LOG2_COPY_RATIO',edgecolor='black',color='dodgerblue',linewidth=0.001,data=df_subset[df_subset['CONTIG']==group],s=0.3,alpha=0.7)
            ax.set_xlabel(group,fontsize=24,rotation=40,ha="right",color='dodgerblue',fontweight='bold')
        elif (group =='chrM'):
            ax.scatter(x='START',y='LOG2_COPY_RATIO',edgecolor='black',color='dodgerblue',linewidth=0.001,data=df_subset[df_subset['CONTIG']==group],s=1,alpha=0.9)
            ax.set_xlabel(group,fontsize=24,rotation=40,ha="right",color='dodgerblue',fontweight='bold')
        else:
            ax.scatter(x='START',y='LOG2_COPY_RATIO',edgecolor='black',color='grey',linewidth=0.001,data=df_subset.loc[df_subset['CONTIG']==group],s=0.3,alpha=0.5)
            ax.set_xlabel(group,fontsize=18,rotation=40,ha="right",color='black',fontweight='bold')
        for index, row in segments.loc[segments['CONTIG']==group].iterrows():
            if row['MEAN_LOG2_COPY_RATIO']>0.33 or row['MEAN_LOG2_COPY_RATIO']<-0.33:
                x,y = (row['START'],row['END']),(row['MEAN_LOG2_COPY_RATIO'],row['MEAN_LOG2_COPY_RATIO'])
                ax.plot(x,y,color='red',linewidth=1)
        sns.despine(left=True)
        ax.set(ylabel=None,facecolor=c_back)
        ax.set_yticks([-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1])
        ax.patch.set_alpha(0.25)
        ax.set_xticks([])
    #ax.axhline(y=row['MEAN_LOG2_COPY_RATIO'], xmin = row['START'], xmax = row['END'])
        sns.despine(left=True)#sns.stripplot(x='Compound',y='value',hue='dataset',edgecolor='black',linewidth=1,data=concat_subset_c,dodge=True,palette=color_pal,size=6,jitter=0.3,ax=ax,alpha=0.9,hue_order = ['Sensitive(NH4)','Sensitive','Resistant'])
        ax.set(ylabel=None,facecolor=c_back)
        ax.patch.set_alpha(0.25)
        ax.set_xticks([])
    #plt.setp(axes, ylim=(-1.5,1.5),yticks=[])
    plt.setp(axes, ylim=(-2,2),yticks=[])
    axes[0].set_yticks([-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5])
    #axes[0].set_yticklabels(['-1.25','-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1','1.25'],fontsize=16,)
    axes[0].set_yticklabels(['-1.5','','-1','','-0.5','','0','','0.5','','1','','1.5'],fontsize=16,)
    all_axes = fig.get_axes()
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
    #        #ax.spines['left'].set_position(('outward',10))
        if ax.is_last_col():
            ax.spines['right'].set_visible(True)
    fig.text(0.035, 0.5, r'Log$_{2}$ Copy Ratio',fontsize=18,fontname="Arial",fontweight='bold',ha='center', va='center', rotation='vertical',color='black',wrap=True)
    plt.savefig('CNV_'+frame_name+'_023123_Final_A.png',dpi=600,transparent=False, bbox_inches='tight', pad_inches=1)

for frame,seg,name in zip(frames,segs,names):
    plot_CNV(frame,seg,name)