#Plot Analyzed CTG
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
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics import jaccard_similarity_score
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
import lmfit
plt.rcParams['font.family']='Arial'
plt.rc('axes', linewidth=2)
import os
import pandas as pd
import sys
import glob
import fnmatch
from itertools import cycle
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14,'axes.labelweight':'bold','legend.handlelength': 1})

#frame = pd.read_csv('Combined_Frame_121121.csv')
#lis1 =['Paclitaxel','B508']
#frame = frame.loc[frame['compound'].isin(lis1)]
#evaluated = pd.read_csv('Combined_Evaluated_Frame_121121.csv')
#print(evaluated.head())
#expanded frames
frame = pd.read_csv('Combined_Frame_063022.csv')
lis1 =['Paclitaxel','B508']
frame = frame.loc[frame['compound'].isin(lis1)]
evaluated = pd.read_csv('Combined_Evaluated_Frame_063022.csv')
#print(evaluated.head())


compoundData = frame.groupby(['compound','frame_name'],sort=False)
#colors = ["steelblue","crimson", "gold"]
colors = ["lightskyblue","tomato"]
alt_colors = ['dodgerblue','crimson']
#cmap = matplotlib.cm.get_cmap('RdYlBu')
cmap = matplotlib.cm.get_cmap('RdBu')
#color_pal = sns.color_palette('RdYlBu', n_colors=3)
#color_pal = sns.color_palette(colors, n_colors=3)
color_pal = sns.color_palette(colors, n_colors=2)
alt_color_pal = sns.color_palette(alt_colors, n_colors=(len(compoundData['frame_name'].unique())-1))
hex_colors = list(color_pal.as_hex())
alt_hex_colors = list(alt_color_pal.as_hex())
nums=[0,1]
licycle = cycle(nums)

#fig,ax = plt.subplots(figsize = (5.5,4.5))
fig,ax = plt.subplots(figsize = (4.75,4))

for name,group in compoundData:
	h=next(licycle)
	num = 0
	if name[1] == 'TRno1.CSV':
		new_label = 'Parental '+str(name[0])
		#ax.scatter(x='Log Concentration',y='Normalized',data=group,color=hex_colors[h],cmap='RdYlBu',s=1)
		ax.plot(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],evaluated['new'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],alpha=1,lw=4,color=alt_hex_colors[h],label=new_label)
		#sns.lineplot(x='dose',y='new',data=evaluated_frame,hue='compound',palette=color_pal)
		#print(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])])
		#plt.fill_between(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],evaluated['lowerbound'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],evaluated['upperbound'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],alpha=0.01,color=alt_hex_colors[h],cmap='RdBu')
	else:
		num = num+1
		new_label = 'Resist Clone '+str(name[1])+' '+str(name[0])
		#ax.scatter(x='Log Concentration',y='Normalized',data=group,color=hex_colors[h],cmap='RdYlBu',s=1)
		ax.plot(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],evaluated['new'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],alpha=0.5,color=hex_colors[h],label=new_label)
		#sns.lineplot(x='dose',y='new',data=evaluated_frame,hue='compound',palette=color_pal)
		#print(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])])
		#plt.fill_between(evaluated['dose'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],evaluated['lowerbound'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],evaluated['upperbound'].loc[(evaluated['compound']==name[0])&(evaluated['frame_name']==name[1])],alpha=0.01,color=hex_colors[h],cmap='RdBu')
plt.xlabel('Log[I], M',fontsize=18,fontname="Arial",fontweight='bold')
plt.xticks(fontsize = 14)
plt.ylabel('Relative Viability',fontsize=18,fontname="Arial",fontweight='bold')
plt.yticks(fontsize = 14)
#plt.ylim(0.15,1.6)
plt.ylim(0.0,1.5)
#plt.xlim(-9.8,-5)
plt.xlim(-9.8,-5)

handles, labels = plt.gca().get_legend_handles_labels()
handles_new = [handles[-4],handles[-3],handles[-2],handles[-1]]
labels_new = ['Resist Clone 1-40 +B508','Resist Clone 1-40 +Paclitaxel','Parental +B508','Parental +Paclitaxel',]
plt.legend(handles_new, labels_new, ncol=2,bbox_to_anchor=(1.2,0.8),
           borderaxespad=0.)
plt.tight_layout()
plt.savefig('CTG_Output_050723_A.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()