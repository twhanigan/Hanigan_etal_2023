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
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14,'axes.labelweight':'bold','legend.handlelength': 1})

import os
import pandas as pd
import sys
import glob
import fnmatch
from itertools import cycle

#frame = pd.read_csv('Combined_Frame_121121.csv')
#evaluated = pd.read_csv('Combined_Evaluated_Frame_121121.csv')
#compoundData = frame.groupby(['compound','frame_name'],sort=False)
#print(evaluated.head())
frame = pd.read_csv('Combined_Frame_063022.csv')
lis1 =['Paclitaxel','B508']
evaluated = pd.read_csv('Combined_Evaluated_Frame_063022.csv')
print(evaluated.head())
compoundData = frame.groupby(['compound','frame_name'],sort=False)
colors = ["steelblue","crimson", "gold"]
alt_colors = ['dodgerblue','crimson','orange']
cmap = matplotlib.cm.get_cmap('RdYlBu')
color_pal = sns.color_palette(colors, n_colors=3)
alt_color_pal = sns.color_palette(alt_colors, n_colors=(len(compoundData['frame_name'].unique())))
hex_colors = list(color_pal.as_hex())
alt_hex_colors = list(alt_color_pal.as_hex())
hex_colors = list(color_pal.as_hex())

nums=[0,1,2]
licycle = cycle(nums)

#fig,ax = plt.subplots(figsize = (5.5,4.5))
fig,ax = plt.subplots(figsize = (4.5,3.75))
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
plt.ylim(0.15,1.6)
#plt.xlim(-10,-5.5)
plt.xlim(-10,-5)

handles, labels = plt.gca().get_legend_handles_labels()
handles_new = [handles[0],handles[1],handles[2],handles[3],handles[4],handles[5]]
labels_new = ['Resist Clone 1-40 +B508','Resist Clone 1-40 +Paclitaxel','Resist Clone 1-40 +5-FU','Parental +B508','Parental +Paclitaxel','Parental +5-FU']
plt.legend(handles_new, labels_new, ncol=2,bbox_to_anchor=(1.2,0.8),
           borderaxespad=0.)
plt.savefig('CTG_Output_063022_5-FU_expanded_C.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
plt.show()