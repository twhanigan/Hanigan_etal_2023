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
import lmfit
from matplotlib import style
from matplotlib.colors import ListedColormap
plt.rcParams['font.family']='Arial'
plt.rc('axes', linewidth=2)

Total=pd.read_csv('Image.csv')
grouped = Total.groupby('Group_Index')



def func(x, a, b, c):
    return a * np.exp(-b * x) + c

#Make Wells
wells = ['%s%02d' % (r, c) for r in 'ABCDEF' for c in range(1, 13)]

#Get plate Layout for Compounds
Neg_Control = ['A01','B01','C01','D01','E01','F01','G01','A12','B12','C12','D12','E12','F12','G12']
Cpd1_0 = ['B01','C01','B12','C12']
Cpd1_1 = ['B02','C02']
Cpd1_2 = ['B03','C03']
Cpd1_3 = ['B04','C04']
Cpd1_4 = ['B05','C05']
Cpd1_5 = ['B06','C06']
Cpd1_6 = ['B07','C07']
Cpd1_7 = ['B08','C08']
Cpd1_8 = ['B09','C09']
Cpd1_9 = ['B10','C10']
Cpd1_10 = ['B11','C11']
Cpd1 = [Cpd1_1,Cpd1_2,Cpd1_3,Cpd1_4,Cpd1_5,Cpd1_6,Cpd1_7,Cpd1_8,Cpd1_9,Cpd1_10,Cpd1_0]

Neg_Control = ['A01','B01','C01','D01','E01','F01','G01','A12','B12','C12','D12','E12','F12','G12']
Cpd2_0 = ['D01','E01','D12','E12']
Cpd2_1 = ['E02','D02']
Cpd2_2 = ['E03','D03']
Cpd2_3 = ['E04','D04']
Cpd2_4 = ['E05','D05']
Cpd2_5 = ['E06','D06']
Cpd2_6 = ['E07','D07']
Cpd2_7 = ['E08','D08']
Cpd2_8 = ['E09','D09']
Cpd2_9 = ['E10','D10']
Cpd2_10 = ['E11','D11']
Cpd2 = [Cpd2_1,Cpd2_2,Cpd2_3,Cpd2_4,Cpd2_5,Cpd2_6,Cpd2_7,Cpd2_8,Cpd2_9,Cpd2_10,Cpd2_0]

Cpd3_0 = ['G01','F01','G12','F12']
Cpd3_1 = ['G02','F02']
Cpd3_2 = ['G03','F03']
Cpd3_3 = ['G04','F04']
Cpd3_4 = ['G05','F05']
Cpd3_5 = ['G06','F06']
Cpd3_6 = ['G07','F07']
Cpd3_7 = ['G08','F08']
Cpd3_8 = ['G09','F09']
Cpd3_9 = ['G10','F10']
Cpd3_10 = ['G11','F11']
Cpd3 = [Cpd3_1,Cpd3_2,Cpd3_3,Cpd3_4,Cpd3_5,Cpd3_6,Cpd3_7,Cpd3_8,Cpd3_9,Cpd3_10,Cpd3_0]

Cpd4_0 = ['G01','H01']
Cpd4_1 = ['G02','H02']
Cpd4_2 = ['G03','H03']
Cpd4_3 = ['G04','H04']
Cpd4_4 = ['G05','H05']
Cpd4_5 = ['G06','H06']
Cpd4_6 = ['G07','H07']
Cpd4_7 = ['G08','H08']
Cpd4_8 = ['G09','H09']
Cpd4_9 = ['G10','H10']
Cpd4_10 = ['G11','H11']
Cpd4 = [Cpd4_1,Cpd4_2,Cpd4_3,Cpd4_4,Cpd4_5,Cpd4_6,Cpd4_7,Cpd4_8,Cpd4_9,Cpd4_10,Cpd4_0]




Control_total = Total.loc[Total['Metadata_Well'].isin(Neg_Control)]
c_0 =[]
c_1 = []
c_2 = []
c_3 = []
c_4 = []
c_5 = []
c_6 = []
c_7 = []
c_8 = []
c_9 = []
c_10 = []
empty = [c_0,c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10]

Compounds = [Cpd1,Cpd2,Cpd3]
#print(Total['Count_Cell'].loc[(Total['Group_Index']==2)&(Total['Metadata_Well'].isin(Cpd1_1))])

def ll4(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''
    return(c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))

# Make density plots for each cell line
Max = Total['Count_Cell'].max()
Min = Total['Count_Cell'].min()
Times = list(Total['Metadata_Time'].unique())
print(Times)

Concentrations = ['10 uM B508','3 uM B508','1 uM B508','0.3 uM B508','0.1 uM B508','0.03 uM B508','0.01 uM B508','0.003 uM B508','0.001 uM B508','0.00033 uM B598','Control']
Conc = [0.00001,3.3333e-6,1.11111e-6,3.7037e-7,1.2345e-7,4.1152e-8,1.3717e-8,4.5724e-9,1.5241e-9,5.08052e-10,0.0]
LogConc = np.log10(Conc)
bins = int((Max-Min)/2)


Compound_Names = ['B508', 'Leflunomide', '5-FU']


def make_frame(cpd):
	d=-1
	lst = []
	for p in cpd:
		c = -1
		d=d+1
		name = Compound_Names[d]
		for i in p:
			c = c+1
			label = Concentrations[c]
			for t in Times:
				empty[c]=pd.DataFrame(Total['Count_Cell'].loc[(Total['Metadata_Well'].isin(i))&(Total['Metadata_Time']==t)])
				mean_intensity = empty[c].mean()
				x_value = Conc[c]
				log_x_value = LogConc[c]
				#Normalized = Total['Response'].loc[frame['compound']==name&frame['Time']==t]-frame['Response'].loc[(frame['compound']==name)&(frame['Time']==t0)]
				lst.append([name,x_value,log_x_value,mean_intensity[0],t])
	final = pd.DataFrame(lst, columns = ['compound','Concentration','Log Concentration','Response','Time'])
	return(final)

def sigmoid(x,b,c,d,e):
    '''This function is basically a copy of the LL.4 function from the R drc package with
     - b: hill slope
     - c: min response
     - d: max response
     - e: EC50'''
    return(c+(d-c)/(1+np.exp(b*(np.log(e)-np.log(x)))))
def exponential(x, a, k,b):
    return a*np.exp(x*k)-b
def expon_simple(a,x,k):
    return a*np.exp(x*k)
def IC50(x,a,b,c):
	# a: top
	# b: bottom
	# c: IC50
	return(b+(a-b)/1+(10**(x-np.log(c))))  

def pDose(x):
	return(np.log10(x))

def normalize(frame):
	t0 = min(Times)
	compoundData = frame.groupby(['compound','Concentration'],sort=False)
	for name,group in compoundData:
		val_min=group.loc[group['Time']==0]
		#print(group['Response'])
		Normalized = []
		#Normalized=group['Response']-val_min['Response'].values
		for row in group.iterrows():
			Normalized=group['Response']-val_min['Response'].values
			frame.loc[(frame['compound']==name[0])&(frame['Concentration']==name[1]),'Normalized']=Normalized
	final = frame.copy()
	return(final)


def exponent(frame):
	compoundData = frame.groupby(['compound','Concentration'],sort=False)
	f=-1
	fitData = []
	evaluated_result =[]
	evaluated_frame = pd.DataFrame(columns = ['time','new','compound','concentration','growth_rate'],dtype=np.float)
	for name,group in compoundData:
		gmodel = Model(exponential)
		t0 = group['Normalized'].loc[group['Time']==0].mean()
		tfinal = group['Normalized'].loc[group['Time']==max(Times)].mean()
		init_slope = (tfinal-t0)/10
		#print(init_slope,t0,tfinal)
		params = gmodel.make_params()
		params.add('a', min=2, max=3)
		params.add('k')
		params.add('b',value=0)
		f=f+1
		xdata = group['Time']
		xeval=np.arange(0,max(Times),1)
		x_logdata = group['Log Concentration']
		ydata = group['Normalized']
		#refDose = np.linspace(min(xdata)*0.95,max(xdata)*1.25,10000)
		#dose = pDose(refDose)
		try:
			result = gmodel.fit(ydata,params,x=xdata,k=init_slope,b=0)
			#result = gmodel.fit(ydata,params,a=2,x=xdata,k=init_slope)
			k_out=result.params.valuesdict()
			new = result.eval(x=xeval)
			dely = result.eval_uncertainty(x=xeval)
			upper_bound = new+dely
			lower_bound = new-dely
			evaluated_frame = evaluated_frame.append(pd.DataFrame({'time':xeval,'new':new,'compound':name[0],'concentration':name[1]}),sort=False)
			frame.loc[(frame['compound']==name[0])&(frame['Concentration']==name[1]),'Growth_Rate']=k_out['k']

		except ValueError:
			result = gmodel.fit(ydata,params,x=xdata,k=-init_slope,b=-1)
			#result = gmodel.fit(ydata,params,a=-2,x=xdata,k=init_slope)
			new = result.eval(x=xeval)
			dely = result.eval_uncertainty(x=xeval)
			upper_bound = new+dely
			lower_bound = new-dely
			evaluated_frame = evaluated_frame.append(pd.DataFrame({'time':xeval,'new':new,'compound':name[0],'concentration':name[1]}),sort=False)
			frame.loc[(frame['compound']==name[0])&(frame['Concentration']==name[1]),'Growth_Rate']=k_out['k']
	#frame.loc[(frame['compound']==name[0])&(frame['Concentration']==name[1]),'Normalized']=Normalized
	final = frame.copy()
	evaluated_frame.to_excel('Evaluated_Growth_Rates.xlsx')
	#print(final)
	return(final)

def Norm_Control(frame):
	compoundData = frame.groupby(['compound'],sort=False)
	for name,group in compoundData:
		val_min = group['Growth_Rate'].loc[group['Concentration']==0].mean()
		print("val_min",val_min)
		Normalized = []
		Normalized = (2**(group['Growth_Rate']/val_min))-1
		frame.loc[(frame['compound']==name),'Normalized_Rate']=Normalized
	final = frame.copy()
	print(final)
	final.to_excel('Final_Frame_042423.xlsx')
	return(final)

def get_curve_lm_b(frame):
	compoundData = frame.groupby(['compound','Time'],sort=False)
	colors = ["limegreen","firebrick", "gold"]
	#color_pal = sns.color_palette(colors, n_colors=len(compoundData))
	color_pal = sns.color_palette(colors, n_colors=3)
	fitData = []
	evaluated_result =[]
	evaluated_frame2 = pd.DataFrame(columns = ['dose','new','dely','upperbound','lowerbound','compound','Time'],dtype=np.float)
	f=-1
	for name,group in compoundData:
		gmodel = Model(sigmoid)
		val_max = group['Normalized_Rate'].loc[group['Concentration']==group['Concentration'].min()].mean()
		val_min = group['Normalized_Rate'].loc[group['Concentration']==group['Concentration'].max()].mean()
		init_slope = (val_max-val_min)/10
		print(val_max,init_slope)
		params = gmodel.make_params()
		params.add('b', max =-4,min=-1)
		params.add('d',max=group['Normalized_Rate'].max()*2)
		params.add('c',min=group['Normalized_Rate'].min())
		params.add('e',min=1e-10)
		f=f+1
		xdata = group['Concentration'].unique()
		x_logdata = group['Log Concentration']
		ydata = group['Normalized_Rate']
		refDose = np.linspace(min(xdata)*0.95,max(xdata)*1.25,10000)
		dose = pDose(refDose)
		try:
			result = gmodel.fit(ydata,params,x=xdata,b=0.08,c=0,d=1,e=1e-5)
			new = result.eval(x=refDose)
			dely = result.eval_uncertainty(sigma=3)
			#dely=0.3*val_min
			upper_bound = new+dely
			lower_bound = new-dely
			evaluated_frame2 = evaluated_frame2.append(pd.DataFrame({'dose':dose,'new':new,'dely':dely,'upperbound':upper_bound,'lowerbound':lower_bound,'compound':name[0],'Time':name[1]}),sort=False)
		except ValueError:
			result = gmodel.fit(ydata,params,x=xdata,b=0.08,c=0.5,d=1,e=1e-7)
			new = result.eval(x=refDose)
			#dely=result.eval_uncertainty(sigma=3)
			dely=0.2*val_max
			upper_bound = new+dely
			lower_bound = new-dely
			evaluated_frame2 = evaluated_frame2.append(pd.DataFrame({'dose':dose,'new':new,'dely':dely,'upperbound':upper_bound,'lowerbound':lower_bound,'compound':name[0],'Time':name[1]}),sort=False)
		except RuntimeError:
			result = gmodel.fit(ydata,params,x=xdata,b=-0.08,c=1,d=1,e=1e-7)
			new = result.eval(x=refDose)
			#dely=result.eval_uncertainty(sigma=3)
			dely=0.2*val_min
			upper_bound = new+dely
			lower_bound = new-dely
			evaluated_frame2 = evaluated_frame2.append(pd.DataFrame({'dose':dose,'new':new,'dely':dely,'upperbound':upper_bound,'lowerbound':lower_bound,'compound':name[0],'Time':name[1]}),sort=False)
	fig, ax = plt.subplots(constrained_layout=False,figsize=(5,4))
	sns.scatterplot(x='Log Concentration',y='Normalized_Rate',data=frame.loc[frame['Time']==0],edgecolor="none",sizes=10,lw=2,hue='compound',palette=color_pal,ax=ax)
	#ax.plot(x='dose',y='new',data=evaluated_frame2,hue='compound',palette=color_pal,legend=False)
	colors = dict(zip(evaluated_frame2['compound'].unique(),colors))
	print(colors)
	for name,group in evaluated_frame2.loc[evaluated_frame2['Time']==0].groupby('compound'):
		ax.fill_between(group['dose'],group['lowerbound'],group['upperbound'],alpha=0.1,color=colors[name],cmap=ListedColormap(color_pal.as_hex()))
		ax.plot(group['dose'],group['new'],color=colors[name])
	#plt.ylim(-0.3,1.6)
	handles, labels = ax.get_legend_handles_labels()
	plt.xlabel('Log[I], M',fontsize=18,fontname="Arial",labelpad=2,fontweight='bold')
	plt.xticks(fontsize = 12)
	plt.ylabel('Growth Rate k(c,t)',fontsize=18,labelpad=2,fontname="Arial",fontweight='bold')
	#plt.yticks([-0.3,0.4,1.1],fontsize = 12)
	#ax.set_yticklabels(['-0.5','0','1.0'])
	plt.legend(handles,labels,loc='upper right',bbox_to_anchor=(1.1, 1.1, 1.1, 1.1))
	plt.tight_layout()
	plt.savefig('Growth_Rate_042423b.png',dpi=600,transparent=True,bbox_inches='tight',pad_inches=2)
	plt.show()
	return(frame,evaluated_frame2)

#a,b= exponent(normalize(make_frame(Compounds)))
#a.to_excel('Testing.xlsx')
#b.to_excel('evaluated_frame.xlsx')
b = pd.read_excel('Final_Frame_042223.xlsx')
#c,d = get_curve_lm_b(Norm_Control(exponent(normalize(make_frame(Compounds)))))
c,d = get_curve_lm_b(b)
d.to_excel('evaluated_frame_042423.xlsx')