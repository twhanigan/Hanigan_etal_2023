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
import conorm
#import rpy2
#print(rpy2.__version__)
import os
#os.environ["R_HOME"] = r"C:\Program Files\R\R-4.2.3"
#os.environ["PATH"]   = r"C:\Program Files\R\R-4.2.3\bin\x64" + ";" + os.environ["PATH"]
#import rpy2.robjects
#from rpy2 import numpy2ri, pandas2ri, Formula
#import rpy2.robjects as robjects
#from rpy2.robjects import numpy2ri, pandas2ri, Formula
#from rpy2.robjects.packages import importr
#pandas2ri.activate()
#numpy2ri.activate()

# import R libraries
#DESeq2 = importr('DESeq2')
#edgeR = importr('edgeR')
#Limma = importr('limma')
#stats = importr('stats')

#df = pd.read_csv('Log2_Frame_MedianNorm_TMT.csv')
df = pd.read_excel('proList_2pep.xlsx')
df = df.drop_duplicates(subset='description')
#Log2_Frame = np.log2(
 #   df.loc[:, 'avg int m/z_126.127726':'avg int m/z_131.13818']).copy()
Log2_Frame = df.loc[:, 'avg int m/z_126.127726':'avg int m/z_131.13818'].copy()
Log2_Frame.replace([np.inf, -np.inf], np.nan, inplace=True)
Log2_Frame.index = df['description'].str.split(" ", 1).str.get(0)
indexes = Log2_Frame.index
print(len(indexes))
dup = indexes.drop_duplicates()
print(len(dup))
Log2_Frame = Log2_Frame[~Log2_Frame.index.duplicated(keep='first')]
print(len(Log2_Frame))
Log2_Frame = Log2_Frame.dropna(how='any', axis=0)

Log2_Frame.columns = ['PC3_Rep1','PC3_Rep2', 'PC3_508_Rep1', 'PC3_508_Rep2', 'PC3_ZFHX3_508_Rep1','PC3_ZFHX3_508_Rep2', 'H460_Rep1', 'H460_Rep2', 'H460_508_Rep1', 'H460_508_Rep2']
#Log2_Frame.columns = ['PC3','PC3', 'PC3_508', 'PC3_508', 'PC3_ZFHX3_508','PC3_ZFHX3_508', 'H460', 'H460', 'H460_508', 'H460_508']
Log2_Frame.to_csv('NoNorm.csv')
nf  = conorm.tmm_norm_factors(Log2_Frame)
df_tmm = conorm.tmm(Log2_Frame)
df_tmm_cpm = conorm.cpm(Log2_Frame, norm_factors=nf)
df_tmm_cpm.to_csv('TMT_Proteomics_TMM_040523c.csv')

print(df_tmm_cpm.head())
H460_Frame = df_tmm_cpm[[ 'H460_Rep1', 'H460_Rep2', 'H460_508_Rep1', 'H460_508_Rep2']]
print(H460_Frame.columns)

design_matrix = pd.DataFrame({'condition':['treated'] * 2 + ['untreated'] * 2})
design_matrix.index = H460_Frame.columns
print(design_matrix)
#DE = DE_rpy2(count_matrix=Log2_Frame, design_matrix=design_matrix)
#DE.deseq2()