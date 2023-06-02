# Protein level abundance associated with B819 Sensitivity
import numpy as np
import pandas as pd
import math
import numpy as d
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
from matplotlib.animation import FFMpegWriter
from matplotlib import style
import statistics

#Raw_Frame_Mutations = pd.read_excel('Mutations_4-18-20.xlsx')
Raw_Frame_Mutations = pd.read_excel('G:\\Python\\Mutation Correlation\\Protein Level\\NSCLC_Proteomics_Final_B.xlsx')
Raw_Frame_Mutations = Raw_Frame_Mutations.drop_duplicates('Site')
Raw_Frame_Mutations.index = Raw_Frame_Mutations['Site']
Raw_Frame_Mutations = Raw_Frame_Mutations.drop('Site',axis=1)
Raw_Frame_Mutations = Raw_Frame_Mutations.drop(['Accession','peptide num','Description'],axis=1).dropna()


#Read in Activity Matrix
sensitivity = pd.read_csv('Activity_Matrix_B508_D.csv')
Activity_data = sensitivity.iloc[0,:].values.tolist()
sens_lines = ['H460' ,'H1975' ,'SW620' , 'H2122', 'A549']
res_lines = ['KP4','MDAMB231','H2009','HCT116','PC3']
Correlation = pd.DataFrame.corrwith(Raw_Frame_Mutations, sensitivity, axis=1)
Correlation.to_csv('Correlation_run2_B.csv')
fig,ax = plt.subplots(figsize=(3.5,5))
Mutation_data = []
Corr_Frame = []
Line_Regress = []
for index,row in Raw_Frame_Mutations.iterrows():
  mut_data = row.iloc[0:].values.tolist()
  Activity_data = sensitivity.iloc[0,:].values.tolist()
  diff  = row.loc[sens_lines].mean() - row.loc[res_lines].mean()
  #print(mut_data,len(mut_data),Activity_data,len(Activity_data))
  slope, intercept, r_value, p_value, std_err = stats.linregress(mut_data,Activity_data)
  Line_Regress.append([index,slope,intercept,r_value,p_value,std_err,diff])
line = pd.DataFrame(Line_Regress, columns = ['index','slope','intercept','r_value','p_value','std_err','diff'])
line.index = line['index']
line = line.drop('index',axis=1)
print(line.head())
final = pd.concat([Raw_Frame_Mutations,line],axis=1)
final.to_csv('Line_Regress_050923.csv')

