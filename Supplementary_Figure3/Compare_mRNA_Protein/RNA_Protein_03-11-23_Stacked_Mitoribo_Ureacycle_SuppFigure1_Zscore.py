import numpy as np
import pandas as pd
import math 
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
from more_itertools import unique_everseen
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
from scipy.optimize import curve_fit
from lmfit.models import LinearModel, StepModel
from lmfit import Model
from matplotlib import gridspec
import re
from adjustText import adjust_text
from pylab import *
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14,'axes.labelweight':'bold'})

#Import and Noramlize All Samples to there protein Concentration
sns.set_style(style="white")
HI = pd.read_excel('proList_MedianNorm_Zscore_1-20-21.xlsx').dropna()
HI.index = HI['Name']
RNA = pd.read_csv('RNA_Log2_Annotated_B.csv').dropna()
RNA.index =RNA['Name']
RNA.loc[:,'KP4':] = stats.zscore(RNA.loc[:,'KP4':],axis=1)
RNA['Name'].loc[RNA['Name']=='COX2'] = 'MT-CO2'
RNA['Name'].loc[RNA['Name']=='SLC25A15'] = 'SLC25A55'
RNA['Name'].loc[RNA['Name']=='SLC2A1'] = 'SLC25A15'
RNA['Name'].loc[RNA['Name']=='CPS1'] = 'EIF4BB'
RNA['Name'].loc[RNA['Name']=='EIF4B'] = 'CPS1'
#RNA['Name'].loc[RNA['Name']=='EIF4B'] = 'CPS1'
HI['Name'].loc[HI['Name']=='PDLIM1'] = 'SLC25A36'
HI['Name'].loc[HI['Name']=='SLC2A1'] = 'SLC25A15'
HI['Name'].loc[HI['Name']=='RAB8A'] = 'NAGS'

mtribo = ['MRPL15','MRPL4','MRPL39','OXA1L']
#ComplexV = ['ATP5G3','ATP5F1','ATP5I','ATP5D','ATP5J2','ATP5O','ATP5L','ATP5B','ATP5H','ATP5C1','ATP5A1','ATP6V1H','MT-ATP6','ATP5I','ATP5H','ATP2A2']
#Glycolysis = ['SLC2A1','SLC2A2','SLC2A3','SLC2A4','SLC2A5','HK1','HK2','HK3','GCK','G6PC','GPI','PFKM','PFKL','PFKP','TPI1','GAPDH','PGAM1','ENO1','ENO2','ENO3','PKM2','PKLR','LDHA','LDHB','LDHC','MPC1','MPC2','GOT1','GOT2','MDK2']
urea = ['ASS1','PYCR2','PYCR1','ALDH18A1','ATP1B3','CPS1','SLC25A15','OAT','NAGS']
pyrim = ['DTYMK','TK2','TK1','CDA','GDA','NTPCR','CAPN2','SLC25A36','DHODH']
#HI = pd.read_excel('proList_MedianNorm_10-20-20.xlsx').dropna().drop('Unnamed: 0',axis=1)

#Log2 Normalize and subtract control

#comp = ['ASS1','CPS1','ORNT','NAGLU','CAD','ASNS','ALDH18A1','PYCR2','PYCR1','ARG2','ASL','ALDH4A1','GOT1','GOT2','SLC25A13','DHODH','CA2','CKB','UMPS','DTYMK']
#comp = ['ASS1','ALDH4A1','CPS1','PYCR2','UMPS','NAGK','PYCR1','ALDH18A1','CKB','PDK3','BCAT1','GDA','OAT','CAD','GLUD1','ARG1','GOT1','DTYMK']
#large urea + mito ribo
#comp_rna = ['MRPL12','MRPL11','MRPL4','MRPL53','MRPL30','DDX28','TBRG4','OXA1L','MT-CO2','MT-CYB','ASS1','CPS1','CDA','DTYMK','NAGS']
#comp = ['MRPL12','MRPL11','MRPL4','MRPL53','MRPL30','DDX28','TBRG4','OXA1L','MT-CO2','MT-CYB','ASS1','CPS1','CDA','DTYMK','SLC7A7']

comp_rna = ['ASS1','CPS1','ALDH18A1','PYCR2','SLC25A15','MRPL15','MRPL4','MRPL39','OXA1L','GDA','CDA','DTYMK','NAGS','TWF2']
comp = ['ASS1','CPS1','ALDH18A1','PYCR2','SLC25A15','MRPL15','MRPL4','MRPL39','OXA1L','GDA','CDA','DTYMK','NAGS','SLC25A36']

names = ['ASS1','CPS1','ALDH18A1','PYCR2','SLC25A15','MRPL15','MRPL4','MRPL39','OXA1L','GDA','CDA','DTYMK','DHODH','SLC25A36']#NDUFB11



#comp_rna = ['MCU','AGK','PHB','GADD45GIP1','ERAL1','TBRG4','MRPL11','MRPL51','ERAL1','MRPL4','MRPL21','MRPS18A','SSBP1','IMMT','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
#comp = ['MCU','AGK','PHB','GADD45GIP1','ERAL1','TBRG4','MRPL11','MRPL51','ERAL1','MRPL4','MRPL21','MRPS18A','SSBP1','IMMT','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
#comp_rna = ['ERAL1','TBRG4','MRPL11','MRPL51','MRPL4','MRPL21','MRPS18A','SSBP1','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
#comp = ['ERAL1','TBRG4','MRPL11','MRPL51','MRPL4','MRPL21','MRPS18A','SSBP1','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
#comp_rna = ['MRPL15','MRPL27','MRPL4','ERAL1','TBRG4','MRPL11','MRPL51','MRPL21','MRPS18A','SSBP1','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']


HI_subset = HI[HI['Name'].isin(comp)].copy()
RNA_subset = RNA[RNA['Name'].isin(comp_rna)].copy()

#Concatenate the frames


#Order the columns
sample_order = comp
col_order = ['KP4','MDAMB231','H2009','H2122','PC3','H460','H1975','SW620','HCT116','A549']
dataset_order = ['Protein','RNA']
sens = ['H460','H1975','SW620','HCT116','A549']
res = ['KP4','MDAMB231','H2009','H2122','PC3']
HI_subset = HI_subset[col_order].copy().unstack(level=0).reset_index()
HI_sens = HI_subset[HI_subset['level_0'].isin(sens)].copy()
HI_res = HI_subset[HI_subset['level_0'].isin(res)].copy()
concat_subset_protein = pd.concat([HI_res.assign(dataset='Resistant'),HI_sens.assign(dataset='Sensitive')],join='inner', axis=0).dropna()
concat_subset_protein.rename(columns={'dataset':'Activity'}, inplace=True)
concat_subset_protein.rename(columns={0:'Value'}, inplace=True)
concat_subset_protein = concat_subset_protein[~(concat_subset_protein.Value == 0)]
#concat_subset_protein.index = concat_subset_protein['Nam'e]


RNA_subset = RNA_subset[col_order].copy().unstack(level=0).reset_index()
RNA_sens = RNA_subset[RNA_subset['level_0'].isin(sens)].copy()
RNA_res = RNA_subset[RNA_subset['level_0'].isin(res)].copy()
concat_subset_RNA = pd.concat([RNA_res.assign(dataset='Resistant'),RNA_sens.assign(dataset='Sensitive')],join='inner', axis=0).dropna()
concat_subset_RNA.rename(columns={'dataset':'Activity'}, inplace=True)
concat_subset_RNA.rename(columns={0:'Value'}, inplace=True)
concat_subset_RNA = concat_subset_RNA[~(concat_subset_RNA.Value == 0)]
#concat_subset_c = pd.concat([concat_subset_protein,concat_subset_RNA],keys=['protein','RNA'],names=['Activity'],axis=0).dropna().reset_index()

#This somewhat works
#concat_subset_c = pd.merge(concat_subset_protein.assign(dataset='Protein'),concat_subset_RNA.assign(dataset='RNA'), how='left')
concat_subset_c = pd.concat([concat_subset_protein.set_index('Name').assign(dataset='Protein'),concat_subset_RNA.set_index('Name').assign(dataset='RNA')], join='outer',axis=0)


#concat_subset_c = pd.concat([concat_subset_protein.assign(dataset='protein'),concat_subset_RNA.assign(dataset='RNA')],join='left', axis=0).dropna().reindex()
concat_subset_c['Combined'] = concat_subset_c['Activity'] + ' ' +concat_subset_c['dataset']
concat_subset_c['Name'] = concat_subset_c.index
concat_subset_c['new_name'] = concat_subset_c['Name'] + concat_subset_c['dataset']

def tt_group(df):
    row = pd.DataFrame(columns=df.columns)
    grouped=df.groupby(['dataset',df.index])
    for name,group in grouped:
        Res_group = group[group['Activity']=='Resistant']['Value']
        print(Res_group)
        Sens_group = group[group['Activity']=='Sensitive']['Value']
        t,p = ttest_ind(Res_group,Sens_group)
        index_num =group.index.unique()
        #print(df.loc[(df.index==name[1])&(df['dataset']==name[0])])
        df.loc[(df.index==name[1])&(df['dataset']==name[0]),'pval']=p
    return(df)

tt_group(concat_subset_c)
print(concat_subset_c.head())


sorterIndex = dict(zip(dataset_order, range(len(dataset_order))))
samp_sorterIdex = dict(zip(sample_order,range(len(sample_order))))
concat_subset_c['Tm_Rank'] = concat_subset_c['dataset'].map(sorterIndex)
concat_subset_c['Row_Rank'] = concat_subset_c.index.map(samp_sorterIdex)
concat_subset_c = concat_subset_c.sort_values(by=['Row_Rank'],ascending = True).copy()
concat_subset_c.drop(['Row_Rank','Tm_Rank'], 1, inplace = True)
concat_subset_c= concat_subset_c.rename(columns={'Name': 'name'})
print(concat_subset_c.head())
#concat_subset_c = concat_subset_c.sort_values(by=['pval'],ascending = True).copy()

#colors1 = ["red",'#0652ff']
colors1 = ["red",'steelblue']
colors2= ['red','steelblue']
trans = (0.1, 0.2, 0.5, 0.3)
colors = ["red",'forestgreen','tomato','limegreen']
#colors2 = ['gold','forestgreen']
#colors  ['red','#0652ff','gold','forestgreen']
color_pal1= sns.color_palette(colors1, n_colors=2)
color_pal2 = sns.color_palette(colors2, n_colors=2)
color_pal = sns.color_palette(colors, n_colors=4)
boxprops = {'edgecolor': 'black', 'linewidth': 0.5,'alpha':0.8}
lineprops = {'color': 'black', 'linewidth': 0.5}
medianprops = {'color': 'black', 'linewidth': 0.5}
boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops,
                       'whiskerprops': lineprops, 'capprops': lineprops,'alpha':0.35})
color_pal_list = [colors1,colors2]
#widths = [0.5, 1=0., 1, 1]
width = [1]
length = len(concat_subset_c['name'].unique())
n=-1
heights= [len(concat_subset_c.loc[concat_subset_c['dataset']=='Protein']),len(concat_subset_c.loc[concat_subset_c['dataset']=='RNA'])]
#Plot for protein
spec = dict(height_ratios=heights,width_ratios=width,hspace=0.01,wspace=0.05)
fig,(ax1,ax2) = plt.subplots(ncols=1,nrows=2,sharex=False, sharey=False,figsize=(7,6),gridspec_kw=spec)
#xpos = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
g=sns.violinplot(x='name',y='Value',hue='Activity',data=concat_subset_c.loc[concat_subset_c['dataset']=='Protein'],split=True,saturation=0.5,ax=ax1,palette=color_pal1,linewidth=0.35,dodge=False,width=1,hue_order=['Sensitive','Resistant'],**boxplot_kwargs)
g=sns.stripplot(x='name',y='Value',hue='Activity',edgecolor='black',linewidth=0.5,data=concat_subset_c.loc[concat_subset_c['dataset']=='Protein'],dodge=True,palette=color_pal1,hue_order=['Sensitive','Resistant'],size=8,jitter=0.25,ax=ax1,alpha=0.9)
handles,labels = ax1.get_legend_handles_labels()
sns.despine(left=True)
ax1.axhline(0,linewidth=1, color='black',ls='-')
ax1.axhline(2,linewidth=0.5, color='black',ls='--',alpha=0.5)
ax1.axhline(-2,linewidth=0.5, color='black',ls='--',alpha=0.5)
ax1.get_legend().remove()
ax1.set(xlabel=None,ylabel=None)
#ax1.set_xticklabels('')
#ax1.set_xticklabels(concat_subset_c['name'].loc[concat_subset_c['dataset']=='Protein'].unique(),fontsize=14,ha='right',fontweight='bold',color='black',rotation=45)

#ax1.set_xticklabels(concat_subset_c['name'].loc[concat_subset_c['dataset']=='Protein'].unique(),fontsize=14,ha='right',fontweight='bold',color='black')
kwargs= dict(fontsize=14,fontweight='bold',fontname='Arial')
#ax1.tick_params(axis='x', labelrotation=45,pad=2)
fig.subplots_adjust(bottom=0.25, left=0.2,wspace=0.05)
#ax1.text(0.5,.05,str(concat_subset_c['dataset'].loc[concat_subset_c['dataset']=='Protein']),size=12,ha='center',transform=ax1.transAxes,color='black',alpha=0.8,fontweight='bold')
medians = 4
nobs = round(concat_subset_c.loc[concat_subset_c['dataset']=='Protein'].groupby('name')['pval'].min(),3)
ys = concat_subset_c.loc[concat_subset_c['dataset']=='Protein'].groupby(['name'])['Value'].min()
nobs = [str(x) for x in nobs]
pos = range(len(nobs))
#for tick,label in zip(pos,ax1.get_xticklabels()):
#	ax1.text(pos[tick], (ys[tick]-2), nobs[tick], horizontalalignment='center', fontsize=12 ,color='k')

	#ax.text(0.5,0.05,concat_subset_c['pval'][0],size=12,ha='center',transform=ax.transAxes,color='black',alpha=0.8,fontweight='bold')
ax1.legend(handles[0:2], labels[0:2],bbox_to_anchor=(1, 1, 1, 1), loc='upper right',
          borderaxespad=0.)
#plt.ylim(-3.25,6)

#Plot for RNA
g=sns.violinplot(x='name',y='Value',hue='Activity',data=concat_subset_c.loc[concat_subset_c['dataset']=='RNA'],ax=ax2,split=True,saturation=0.5,palette=color_pal2,linewidth=0.5,dodge=False,width=1,hue_order=['Sensitive','Resistant'],**boxplot_kwargs)
g=sns.stripplot(x='name',y='Value',hue='Activity',edgecolor='black',linewidth=0.5,data=concat_subset_c.loc[concat_subset_c['dataset']=='RNA'],dodge=True,palette=color_pal2,hue_order=['Sensitive','Resistant'],size=8,jitter=0.25,ax=ax2,alpha=0.9)
handles,labels = ax2.get_legend_handles_labels()
sns.despine(left=True)
ax2.axhline(0,linewidth=1, color='black',ls='-')
ax2.axhline(2,linewidth=0.5, color='black',ls='--',alpha=0.5)
ax2.axhline(-2,linewidth=0.5, color='black',ls='--',alpha=0.5)
ax2.get_legend().remove()
ax1.set(xlabel=None,ylabel=None)
ax2.set(xlabel=None,ylabel=None)
#ax2.set_xticks(range(0,8))
#ax2.set_xticklabels(concat_subset_c['name'].loc[concat_subset_c['dataset']=='RNA'].unique(),fontsize=14,ha='right',fontweight='bold',color='black')
ax2.set_xticklabels(names,fontsize=14,ha='right',fontweight='bold',color='black')
ax1.set_xticklabels('')

kwargs= dict(fontsize=14,fontweight='bold',fontname='Arial')
ax2.tick_params(axis='x', labelrotation=45,pad=2)
#fig.subplots_adjust(bottom=0.5, left=0.2,wspace=0.05)
#ax2.text(0.5,.05,str(concat_subset_c['dataset'].loc[concat_subset_c['dataset']=='RNA']),size=12,ha='center',transform=ax2.transAxes,color='black',alpha=0.8,fontweight='bold')
medians = 4
nobs = round(concat_subset_c.loc[concat_subset_c['dataset']=='RNA'].groupby(['name'])['pval'].min(),3)
ys = concat_subset_c.loc[concat_subset_c['dataset']=='RNA'].groupby(['name'])['Value'].min()
nobs = [str(x) for x in nobs.tolist()]
pos = range(len(nobs))
#for tick,label in zip(pos,ax2.get_xticklabels()):
#	ax2.text(pos[tick], (ys[tick]-2), nobs[tick], horizontalalignment='center', fontsize=12, color='k')

	#ax.text(0.5,0.05,concat_subset_c['pval'][0],size=12,ha='center',transform=ax.transAxes,color='black',alpha=0.8,fontweight='bold')
ax2.legend(handles[0:2], labels[0:2],bbox_to_anchor=(1., 1, 1, 1), loc='upper right',
          borderaxespad=0.)
#plt.setp(ax1.collections[::4], alpha=.3)
#plt.setp(ax2.collections[::4], alpha=.3)

all_axes = fig.get_axes()
# show only the outside spines
for ax in all_axes:
    for sp in ax.spines.values():
        sp.set_visible(True)
        #sp.set_visible(False)
#    if ax.is_first_row():
#        #ax.spines['bottom'].set_visible(True)
#        #ax.spines['bottom'].set_position(('outward',10))
#        ax.spines['top'].set_visible(True)
#    if ax.is_last_row():
#        ax.spines['bottom'].set_visible(True)
#        #ax.spines['bottom'].set_position(('outward',10))
#    if ax.is_first_col():
#        ax.spines['left'].set_visible(True)
#        #ax.spines['left'].set_position(('outward',10))
#    if ax.is_last_col():
#       ax.spines['right'].set_visible(True)
ax1.set(xlabel=None,ylabel=None)
colors = ["black",'blue','mistyrose','tomato','red','crimson','darkred','grey']
c=-1
ticks = plt.gca().get_xticklabels()
print(ticks)
for tick in ticks:
    c=c+1
    print(tick)
    if tick.get_text() in urea:
        tick.set_color('crimson')
    elif tick.get_text() in pyrim:
        tick.set_color('black')
    elif tick.get_text() in mtribo:
        tick.set_color('gold')
    else:
        tick.set_color('gray')
plt.tight_layout()
fig.text(0.27, 0.85,'Protein',fontsize=14,fontname="Arial",fontweight='bold',ha='center', va='center', rotation='horizontal',color='black',wrap=True)
fig.text(0.28, 0.53,'Transcript',fontsize=14,fontname="Arial",fontweight='bold',ha='center', va='center', rotation='horizontal',color='black',wrap=True)

fig.text(0.12, 0.55,'Z-Score Log$_{2}$ Expression',fontsize=18,fontname="Arial",fontweight='bold',ha='center', va='center', rotation='vertical',color='black',wrap=True)
plt.savefig('MitoRibo_RNA_Protein_061923_F.png',dpi=300,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()


