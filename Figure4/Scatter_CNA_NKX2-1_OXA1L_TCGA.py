import numpy as np
import pandas as pd
import math
import scipy.stats as stats
import warnings
import seaborn as sns
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

# Set matplotlib style
mplstyle.use('seaborn')

# Set matplotlib parameters
rc = {
    "fontname": "Arial",
    'fontsize': 24,
    'fontweight': 'bold',
    'lines.linewidth': 2
}
plt.rcParams.update({"font.family": 'sans-serif',
                     'font.sans-serif': "Arial",
                     'font.size': 12,
                     'axes.linewidth': 2,
                     'axes.labelsize': 14,
                     'axes.labelweight': 'bold'})

# Read data from a CSV file
df = pd.read_csv('TCGA_NKX2-1_OXA1L_052223.txt', sep='\t', comment='@', skip_blank_lines=True)
# df = pd.read_csv('OXA1L_052223.txt', sep='\t', comment='@', skip_blank_lines=True)
print(df.columns)

# Calculate the logarithm base 2 of a column and assign it to a new column 'log10'
df['log10'] = np.log2(df['OXA1L: mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2) (log2(value + 1))'])

# Create a figure and axes
fig, axes = plt.subplots(figsize=(2, 3.5))

# Define colors for the stripplot
colors = ["steelblue", 'dodgerblue', 'grey', 'red', 'firebrick']
color_pal = sns.color_palette(colors, n_colors=5)

# Create a stripplot
sns.stripplot(x='NKX2-1: Putative copy-number alterations from GISTIC',
              y='log10',
              data=df,
              jitter=0.35,
              edgecolor='black',
              palette=color_pal,
              linewidth=0.1,
              size=3,
              ax=axes)

# Set the y-axis limits and tick labels
plt.ylim(9, 13.5)
plt.yticks(fontsize=14)

# Set the x-axis tick labels
axes.set_xticklabels(['N-2', 'N-1', 'N', 'N+1', 'N+2'],
                     fontsize=14,
                     rotation=40,
                     ha="right",
                     fontweight='bold')

# Set x and y labels
plt.xlabel('NKX2-1 CNA\n(chr14q11)', fontsize=18, fontweight='bold')
plt.ylabel('Relative OXA1L mRNA\n'+r'Log$_{2}$(FPKM+1)',
           fontsize=18, fontweight='bold', labelpad=10)

# Save the figure
plt.savefig('CAMTA1_OXA1L_052223_G.png', dpi=600, transparent=True, bbox_inches='tight', pad_inches=1)

# Show the plot
plt.show()
