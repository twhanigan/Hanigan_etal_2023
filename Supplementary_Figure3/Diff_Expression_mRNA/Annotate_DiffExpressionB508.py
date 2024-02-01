import pandas as pd
from bioinfokit import visuz
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
from sklearn.cluster import KMeans
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
import mygene

Mutations= pd.read_excel('Resistant_Final.xlsx')
mg = mygene.MyGeneInfo()
gene_name_list = Mutations['description']#.str.split("] ", 1).str.get(1).str.split(';',1).str.get(0)
Mutations.index = gene_name_list
protein_name = mg.querymany(gene_name_list, scopes='ensembl.gene', fields=[
                            'symbol','genomic_pos.chr','genomic_pos','summary','exons','kegg','go','pathway','reactome','ensembl.type_of_gene'], species = 'human', as_dataframe = True, df_index = True)
protein_name = protein_name.drop_duplicates(subset='symbol')
merged_mutations = Mutations.merge(protein_name, how='left', left_index=True, right_index=True)
merged_mutations.to_csv('Resistant_Final_Annotated.csv')