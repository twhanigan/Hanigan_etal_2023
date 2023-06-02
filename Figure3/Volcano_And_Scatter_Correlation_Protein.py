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
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14,'axes.labelweight':'bold','legend.handlelength': 1})
plt.rcParams['font.family']='Arial'
plt.rc('axes', linewidth=2)
dfH460 = pd.read_csv('Line_Regress_050923.csv')
Raw_Frame = dfH460[['index','KP4', 'H2009',  'H2122',   'MDAMB231', 'PC3', 'H460',   'H1975',  'SW620',   'HCT116','A549']]
#dfH460 = pd.read_csv('old.csv')
sensitivity = pd.read_csv('Activity_Matrix_B508_D.csv')

dfH460.index = dfH460['index']
dfH460['significant'] = (dfH460['p_value'] < 0.001) & ((dfH460['slope'] > 0.5) | (dfH460['slope'] < -0.5))
subset = dfH460.loc[dfH460['significant']==True]

#sign_line=True,genenames=({'ENSG00000241429':'EEF1A1',"ENSG00000240036":'MT-CO2',"ENSG00000211459":"MT-RNR1","ENSG00000198938":"MT-CO3","ENSG00000114923":"SLC4A3","ENSG00000214026":"MRPL23","ENSG00000159423":"ALDH4A1"}),show=True)

#style.use('seaborn-muted')
Mitorib = ['MCU','AGK','PHB','GADD45GIP1','ERAL1','TBRG4','MRPL11','MRPL51','ERAL1','MRPL4','MRPL21','MRPS18A','SSBP1','IMMT','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
OXA1L = ['OXA1L']
ComplexIV = ['MT-CO1','MT-CO2','MT-CO3','COX5B','COX5A','COX6C','COX7C','MT-CYB','COX6B1','COX7A2','PET100']
ComplexI = ['NDUFC1','NDUFB8','NDUFB6','NDUFB3','NDUFA1','NDUFA12','NDUFS5','NDUFA2','NDUFA8','NDUFA6','NDUFAB1','NDUFV2','NDUFA11','NDUFS4','NDUFS6','NDUFS5','NDUFA5','NDUFB10','NDUFA10']
ComplexII = ['SDHD','SDHA','SDHC','SDHAF2','FH','NDUFV1','UQCRFS1','NDUFS2','SUCLG1','SUCLA2','SDF2']
ComplexIII = ['UQCRB','MT-CYB','UQCRC1','CYC1','UQCRH','NDUFV1','UQCRC2','UQCR10','UQCRFS1','UQCRQ','UQCR11']
ComplexV = ['ATP5G3','ATP5F1','ATP5I','ATP5D','ATP5J2','ATP5O','ATP5L','ATP5B','ATP5H','ATP5C1','ATP5A1','ATP6V1H','MT-ATP6','ATP5I','ATP5H','ATP2A2']
#Glycolysis = ['SLC2A1','SLC2A2','SLC2A3','SLC2A4','SLC2A5','HK1','HK2','HK3','GCK','G6PC','GPI','PFKM','PFKL','PFKP','TPI1','GAPDH','PGAM1','ENO1','ENO2','ENO3','PKM2','PKLR','LDHA','LDHB','LDHC','MPC1','MPC2','GOT1','GOT2','MDK2']
UreaCycle = ['ASS1','PYCR2','PYCR1','ALDH18A1','ATP1B3','CPS1','SLC25A15','OAT']
Thymidine = ['DTYMK','TK2','TK1','CDA','NTPCR','CAPN2','LDHA']


index_ComplexI = dfH460.index.isin(ComplexI)
index_ComplexII = dfH460.index.isin(ComplexII)
index_ComplexIII = dfH460.index.isin(ComplexIII)
index_ComplexIV = dfH460.index.isin(ComplexIV)
index_ComplexV = dfH460.index.isin(ComplexV)
index_UreaCycle = dfH460.index.isin(UreaCycle)
index_Thymidine = dfH460.index.isin(Thymidine)
index_mitorib = dfH460.index.isin(Mitorib)
index_OXA1L = dfH460.index.isin(OXA1L)

#index_Glycolysis = dfH460.index.isin(Glycolysis)

dfH460['Type'] = 'Ather'

dfH460['Type'].loc[index_ComplexI] ='Complex I'
dfH460['Type'].loc[index_ComplexII] ='Complex II'
dfH460['Type'].loc[index_ComplexIII] ='Complex III'
dfH460['Type'].loc[index_ComplexIV] ='Complex IV'
dfH460['Type'].loc[index_ComplexV] = 'Complex V'
#dfH460['Type'].loc[index_Glycolysis] = 'Glycolysis'
#dfH460['Type'].loc[index_UreaCycle] = 'Urea Cycle'
#dfH460['Type'].loc[index_Thymidine] = 'De Novo Thymidine Synthesis'
dfH460['Type'].loc[index_mitorib] ='Mitoribosome'
dfH460['Type'].loc[index_OXA1L] ='OXA1L'


dfH460['Type']=dfH460['Type'].fillna(0)

print(dfH460['Type'])
#sns.set(style='white')

#colors = ["red",'forestgreen','gold','#0652ff','grey']
#colors = ["red",'#0652ff','forestgreen','steelblue','orange','violet','tomato','gold','grey']
colors = ['dodgerblue',"gold",'mistyrose','tomato','red','crimson','darkred','grey']

#color_pal= sns.color_palette(colors, n_colors=5)
#color_pal= sns.color_palette(colors, n_colors=9)
color_pal= sns.color_palette(colors, n_colors=8)

length = len(dfH460['Type'].unique())
print(length)
#custom = sns.palplot(sns.xkcd_palette(colors),3)
#orders = ['Mitoribosome','Urea Cycle','De Novo Thymidine Synthesis','Complex I','Complex II','Complex III','Complex IV','Complex V','Other']
orders = ['OXA1L','Mitoribosome','Complex I','Complex II','Complex III','Complex IV','Complex V','Ather']
#sizes = ['Ather','Mitoribosome','Complex I','Complex II','Complex III','Complex IV','Complex V','De Novo Thymidine Synthesis','Urea Cycle']
sizes = ['Ather','OXA1L','Mitoribosome','Complex I','Complex II','Complex III','Complex IV','Complex V']

dfH460=dfH460.sort_values(by=['Type','p_value',],ascending=[True,False])
#fig,ax = plt.subplots(figsize=(4,6))
#fig,ax = plt.subplots(figsize=(3.5,5))
fig,ax = plt.subplots(figsize=(3.5,3.5))

g=sns.scatterplot(x='slope',y='p_value',data= dfH460,hue='Type',size='Type',hue_order = orders,
            sizes=[10,200,150,100,100,100,100,100], size_order=sizes,alpha=0.5, palette=color_pal,linewidth=0.1,edgecolor='black',ax=ax
            )


#plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.yscale('log')
plt.ylim(1,0.000001)
plt.xlim(-2,2)
#plt.axvline(-1,linewidth=1, color='black',ls='--',zorder=-5)
#plt.axvline(1,linewidth=1, color='black',ls='--',zorder=-5)
#plt.axhline(0.01,linewidth=1, color='black',ls='--',zorder=-5)


#Labels
#labels = ["MRPL12","DDX5","NSUN4","MRPL54","ERAL1","MRPL34","MRPL55","MALSU1","MRPL11","MRPL49","MRPL32","MRPL40","MRPL22","MRPL4","MRPL28","MRPL38","MRPL9","MRPL24","MRPL10","MRPL50","MRPL43","MRPL15","MRPL44","MRPS30","MRPL47"]
#labels = ['ERAL1','TBRG4','MRPL11','MRPL51','ERAL1','MRPL4','MRPL21','MRPS18A','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
labels = ['OXA1L']
#keys = ['EEF1A1','MT-CO2','MT-RNR1','MT-CO3','SLC4A3','MRPL23','ALDH4A1']
subset_labels = dfH460.loc[dfH460.index.isin(labels)].copy()
subset_labels['slope']=subset_labels['slope'].astype(float)

#count = -1
#for i in labels:
#	count = count+1
#	x = subset_labels.slope[i]
#	y = subset_labels.p_value[i]
#	plt.text(x,y,labels[count])

texts = subset_labels.sort_values(by=['p_value'],ascending=[False]).reset_index(drop=True)
texts1 = [plt.text(texts.slope[i], texts.p_value[i],
                   'OXA1L', color='black',fontweight='bold', fontsize=14) for i in range(len(texts))]
#texts1 = [plt.text(subsetsens_sub.fold_change[i], subsetsens_sub.pval[i],
#                   subsetsens_sub.Name[i], color='crimson', fontsize=14) for i in range(len(subsetsens_sub))]
#texts2 = [plt.text(subsetsens_sub.fold_change[i], subsetsens_sub.pval[i],
#                   subsetsens_sub.Name[i], color='crimson', fontsize=14) for i in range(len(subsetsens_sub))]

adjust_text(texts1, arrowprops=dict(arrowstyle='->', color='dodgerblue', alpha=1),ha='right',expand_points=(3,1),force_points=2,expand_text=(0, 3))#)#

#leg = g.legend
#leg.set_bbox_to_anchor([1.1, 0.95])  # coordinates of lower left of bounding box
#leg._loc = 2  # if required you can set the loc
plt.xlabel(r'Log$_{2}$'+'FC Protein Expression\n(Sensitive-Resistant)',fontsize=18)
plt.xticks(fontsize = 12)
plt.ylabel('PValue',fontsize=18)
plt.yticks(fontsize = 12)
#plt.tight_layout()
handles, labels = plt.gca().get_legend_handles_labels()
handles_new = handles[0], handles[1], handles[4]
labels_new = ['OXA1L','Mitoribosome','Mutation']
plt.legend(handles_new, labels_new,bbox_to_anchor=(2, 1),borderaxespad=0,ncol=1)
#plt.text(-0.96,0,'Sensitive',color='red',fontsize=16,fontweight='bold',fontname="Arial",wrap=True)
#plt.text(0.04,0,'Resistant',color='blue',fontsize=16,fontweight='bold',fontname="Arial",wrap=True)

plt.savefig('Volcano_ProteinCorr_slope_050923_diff_linregress_B.png',dpi=300,transparent=True,bbox_inches='tight',pad_inches=2)
#plt.show()


# Import Spreadsheets
stacked = Raw_Frame.set_index(['index']).stack().reset_index()
stacked.columns = ['index','Cell Lines','Log_Expression']
sens_stacked = sensitivity.stack().reset_index()
sens_stacked.columns = ['index','Cell Lines','pGC50']
Final_Table = stacked.merge(sens_stacked,on='Cell Lines')
print(Final_Table.head())
Mitoribo = ['MCU','AGK','PHB','GADD45GIP1','ERAL1','TBRG4','MRPL11','MRPL51','ERAL1','MRPL4','MRPL21','MRPS18A','SSBP1','IMMT','MRPL15','MRPL30','MRPS30','MRPL3','MRPL39','MRPL19','MRPL32','MRPL27','OXA1L','MRPL45','MRPL24','ELAC2','MRPL2','FASTKD2','MRPL37','MRP63','DDX5','MRPL28']
OXA1L = ['OXA1L']
ComplexIV = ['MT-CO1','MT-CO2','MT-CO3','COX5B','COX5A','COX6C','COX7C','MT-CYB','COX6B1','COX7A2','PET100']
ComplexI = ['MT-ND1','MT-ND2','MT-ND3','MT-ND4','ND5','NDUFC1','NDUFB8','NDUFB6','NDUFB3','NDUFA1','NDUFA12','NDUFS5','NDUFA2','NDUFA8','NDUFA6','NDUFAB1','NDUFV2','NDUFA11','NDUFS4','NDUFS6','NDUFS5','NDUFA5','NDUFB10','NDUFA10']
ComplexII = ['SDHD','SDHA','SDHC','SDHAF2','FH','NDUFV1','UQCRFS1','NDUFS2','SUCLG1','SUCLA2','SDF2']
ComplexIII = ['UQCRB','MT-CYB','UQCRC1','CYC1','UQCRH','NDUFV1','UQCRC2','UQCR10','UQCRFS1','UQCRQ','UQCR11']
ComplexV = ['MT-ATP6','MT-ATP8','ATP5G3','ATP5F1','ATP5I','ATP5D','ATP5J2','ATP5O','ATP5L','ATP5B','ATP5H','ATP5C1','ATP5A1','ATP6V1H','MT-ATP6','ATP5I','ATP5H','ATP2A2']

colors = ["gold",'dodgerblue',]
#color_pal= sns.color_palette(colors, n_colors=5)
#color_pal= sns.color_palette(colors, n_colors=9)
print(len(Final_Table['index_x'].loc[Final_Table['index_x'].isin(Mitoribo)].unique()))
color_pal= sns.color_palette('blend:yellow,goldenrod', n_colors=31)
#g=sns.lmplot(x='pGC50',y='Log_Expression',hue='index_x',palette=color_pal, data= Final_Table[Final_Table['index_x'].isin(Mitoribo)],facet_kws = dict(sharex=False, sharey=True),scatter_kws = {'alpha':0.7,'linewidth':0.1,'edgecolor':'black'})#'vlag'
fig,ax = plt.subplots(figsize=(3.5,3.5))
sns.regplot(x='pGC50',y='Log_Expression',scatter=False,color='k',data= Final_Table[Final_Table['index_x'].isin(Mitoribo)],line_kws={'alpha':0.6,'lw':2,'zorder':-2},ax=ax)#**kwargs3
sns.scatterplot(x='pGC50',y='Log_Expression',hue='index_x',data= Final_Table[Final_Table['index_x'].isin(Mitoribo)],alpha=0.7,linewidth=0.1,s=50,edgecolor='black',palette=color_pal,ax=ax)#**kwargs3
sns.scatterplot(x='pGC50',y='Log_Expression',data= Final_Table[Final_Table['index_x']=='OXA1L'],alpha=0.7,linewidth=0.1,s=50,edgecolor='black',color='dodgerblue',ax=ax)#**kwargs3

handles, labels = plt.gca().get_legend_handles_labels()
#ax.get_legend().remove()
plt.legend(handles, labels, title='Compound',
           bbox_to_anchor=(1.33, 1.12), fontsize=10)
def annotate(data, **kws):
        m, b, r_value, p_value, std_err = \
            stats.linregress(data['pGC50'],data['Log_Expression'])
        #ax.text(0.05, 1.05, data['Compound'].values[0],fontweight='bold',fontsize='16',color=kws.get("color","k"),transform=ax.transAxes)
        ax.text(0.05, 0.9, f"Slope = {m:.2f}x",fontweight='bold',fontsize='14',
                transform=ax.transAxes)
        ax.text(0.05, .82, f'r^2 ={r_value**2:.3f}',fontweight='bold',fontsize='14',
                transform=ax.transAxes)

annotate(Final_Table[Final_Table['index_x'].isin(Mitoribo)])
#g.set_xticklabels(['-2','-1','0','1','2'])
#g.set_yticklabels(['0','20','40','60','80','100'])
ax.set_ylabel(None)
ax.set_xlabel(None)
ax.set_title('')
ax.tick_params(labelleft=True)
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)
fig.text(0.5,0.01, r'pGC$_{50}$', ha="center", va="center",fontsize=18,fontname="Arial",fontweight='bold')
fig.text(0.0,0.5, r'Log$_{2}$'+'Protein Expression', ha="center", va="center",fontsize=18,fontname="Arial",fontweight='bold',rotation=90)
        #g.fig.subplots_adjust(hspace=0.3,wspace=0.2)
#plt.gca().set_xlim(9,11.5)
#plt.gca().set_xticks([9,10,11,12])
#g.axes[1].set_xlim(9,12)
#g.axes[1].set_xticks([9,10,11,12])
#plt.ylim(0,100)
plt.savefig('Scatterplot_mitoribo_Diff_linregress_A.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()

