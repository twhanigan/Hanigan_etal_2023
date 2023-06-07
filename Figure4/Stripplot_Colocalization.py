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
from scipy.optimize import curve_fit
from lmfit.models import LinearModel, StepModel
from lmfit import Model
from sklearn.utils import shuffle
import lmfit
from matplotlib import style
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':14,'axes.linewidth':2,'axes.labelsize':14,'axes.labelweight':'bold'})



###Import Spreadsheets
Total = pd.read_excel('Aggregated_Zstack_Colocalization_C.xlsx')


colors = ["red",'dodgerblue']
color_pal= sns.color_palette(colors, n_colors=2)

boxprops = {'edgecolor': 'black', 'linewidth': 1}
lineprops = {'color': 'black', 'linewidth': 1}
medianprops = {'color': 'black', 'linewidth': 1}
boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops,
                       'whiskerprops': lineprops, 'capprops': lineprops})
other_kwargs = dict({'edgecolor': 'w', 'linewidth': 1.5,'errcolor':'black','errwidth':2,'capsize':0.15,'line_kws':{'alpha':0.1}})

heights = [2,2]
spec = dict(height_ratios=heights,hspace=-1,wspace=0)
fig,ax = plt.subplots(figsize=(1.1,3))
#g=sns.violinplot(x='B508',y='Intensity_Brdu',hue='B508',data=subsetb[subsetb['B508']=='resistant'],dodge=True,palette='Blues',ax=ax2,inner=None,alpha=0.1)
g=sns.barplot(x='Cell_Line',y='Correlation',hue='Treatment',linewidth=1,data=Total,alpha=0.75,palette=color_pal,ci='sd')#hue_order=labels_new,)
g=sns.stripplot(x='Cell_Line',y='Correlation',hue='Treatment',data=Total,dodge=0.01,palette=color_pal,alpha=0.7,edgecolor='k',linewidth=0.5,size=8,ax=ax,jitter=0.25)
#g=sns.violinplot(x='B508',y='Intensity_Brdu',hue='B508',data=subsetb[subsetb['B508']=='sensitive'],dodge=True,palette='Reds',ax=ax1,inner=None,alpha=0.1)
#g=sns.stripplot(x='Concentration',y='Intensity_Brdu',hue='dataset',data=subsetb[subsetb['B508']=='sensitive'],dodge=True,palette='Reds',alpha=0.45,edgecolor='k',linewidth=0.1,ax=ax3,jitter=0.25)

#g1=sns.pointplot(x='Cell_Line',y='Correlation',hue='Treatment',data=Total,
#              dodge=0.2,
#              join=False, color='k',alpha=0.8,
#              markers="_", scale=1.5, ci=None, ax=ax)


#g=sns.pointplot(x='Time',y='Value',hue='Compound',data=df_new[df_new['Cell_Line']=='H460'],ax= ax,palette=color_pal,sort=False,lw=2,**other_kwargs,legend=False,markers='.')
#g=sns.pointplot(x='Time',y='Value',hue='Compound',data=df_new[df_new['Cell_Line']=='PC3'],ax= ax,palette='Blues',sort=False,lw=2,**other_kwargs,legend=False)
ax.get_legend().remove()

#plt.setp([g.get_children()],alpha=0.4)
#plt.setp([g.lines],alpha=0.5)

#sns.scatterplot(x='Time',y='Average',hue='Compound',edgecolor='black',linewidth=1,data=df.loc[df['Cell_Line']=='H460'],ax=ax,palette='RdBu',alpha=0.9)
#sns.lineplot(x='Time',y='Average',hue='Compound',data=df.loc[df['Cell_Line']=='H460'],ax= ax,palette='RdBu',sort=False,alpha=0.4)
#plt.errorbar(x=df['Time'], y=df['Average'], fmt='none', yerror=df['Stdev'], ecolor='k', elinewidth=2)
#sns.scatterplot(x='Time',y='Normalized',hue='Concentration',edgecolor='black',linewidth=1,data=concat_c,ax=ax,palette='RdBu',alpha=0.9)
#sns.lineplot(x='Time',y='Normalized',hue='Concentration',data=concat_c,ax= ax,palette='RdBu',markers=False,ci=False,sort=False,alpha=0.4)

handles,labels = ax.get_legend_handles_labels()
labels_new = ['Control','+B508']

kwargs= dict(fontsize=14,fontweight='bold',fontname='Arial')

#sns.catplot(x='Sample',y='value',hue='dataset',data=concat_subset_c.loc[concat_subset_c['Compound']==compound],ax= ax, kind="point",dodge=True,palette='RdBu',ci=None)
#sns.lineplot(x='Compound',y='value',hue='dataset',data=concat_subset_c,ax= ax,palette=color_pal,markers=False,ci=False,sort=False,alpha=0.4)
#sns.despine(left=True)
#ax.axhline(0,linewidth=1, color='black',ls='-')
#ax.axhline(2,linewidth=0.5, color='black',ls='--',alpha=0.5)
#ax.axhline(-2,linewidth=0.5, color='black',ls='--',alpha=0.5)
kwargs= dict(fontsize=18,fontweight='bold',fontname='Arial')
plt.xticks(ha="center",fontsize=16)#rotation=45
#ax.set_xticklabels(['S423InsAGC','S423S'],fontsize=16,fontweight='bold',rotation='45',ha='right',fontname='Arial')
ax.set_xticklabels([r'S419(AGC)$_{4}$',r'S419(AGC)$_{5}$'],fontsize=16,fontweight='bold',rotation='45',ha='right',fontname='Arial')
ticks = ax.xaxis.get_majorticklabels()
ticks[1].set_color('red')
ticks[0].set_color('steelblue')

#sns.barplot(x='value',y='Compound',hue='dataset',data=concat_subset_c,palette=color_pal,dodge=0.1,**other_kwargs)
#sns.barplot(x='Compound',y='value',hue='dataset',data=concat_subset_c,palette=color_pal,dodge=True,ax=ax,**other_kwargs)

#g=sns.swarmplot(x='Concentration',y='Intensity_Brdu',hue='Concentration',data=subsetb,palette=color_pal,alpha=0.2,edgecolor=None,linewidth=None,ax=ax)

plt.legend(handles[0:2], labels_new,bbox_to_anchor=(1.5,1.5))
#          borderaxespad=0.)
plt.ylim(0.2,1)
all_axes = fig.get_axes()
# show only the outside spines
for ax in all_axes:
    for sp in ax.spines.values():
        sp.set_visible(False)
        #sp.set_visible(False)
    if ax.is_first_row():
        ax.spines['bottom'].set_visible(True)
#        #ax.spines['bottom'].set_position(('outward',10))
#        ax.spines['top'].set_visible(True)
#    if ax.is_last_row():
#        ax.spines['bottom'].set_visible(True)
#        #ax.spines['bottom'].set_position(('outward',10))
    if ax.is_first_col():
        ax.spines['left'].set_visible(True)
#        #ax.spines['left'].set_position(('outward',10))
#    if ax.is_last_col():
#       ax.spines['right'].set_visible(True)

plt.ylabel('Manders Correlation\n(OXA1L/MT-CO2)',fontsize=18,fontname="Arial",fontweight='bold')

plt.savefig('Colocalization_Strip_071422_Finalh.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()
