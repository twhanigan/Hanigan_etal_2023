import math
import numpy as np
import pandas as pd
from pandas import ExcelWriter, ExcelFile
import scipy
from scipy import stats
from scipy.stats import pearsonr
import warnings
from scipy.stats import ttest_ind
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from adjustText import adjust_text
import random as randy
from matplotlib import style
import statistics
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Set plot style and font properties
rc = {"fontname": "Arial", 'fontsize': 24, 'fontweight': 'bold', 'lines.linewidth': 2}
plt.rcParams.update({"font.family": 'sans-serif', 'font.sans-serif': "Arial", 'font.size': 12,
                     'axes.linewidth': 2, 'axes.labelsize': 14})
plt.rc('axes', linewidth=2)

# Import Spreadsheets
man_df = pd.read_excel('Combined_RawData_081522_percell_b.xlsx').reset_index(drop=True)
Activity = pd.read_csv('Activity_Matrix_B508_D.csv').T.reset_index()
Activity.columns = ['Cell Line', 'pGI50', 'Fold_Change_ATP', 'B']
Activity['GI50'] = 1 * 10 ** -Activity['pGI50']
man_df = man_df.merge(Activity, on='Cell Line', how='left')
man_df = man_df.loc[(man_df['dataset'] == 'Normalized_OCR')]
grouped_df = man_df.groupby(['Cell Line', 'Medium', 'Treatment'], as_index=False)
min_time = man_df['Time'].min()
min_val = man_df.loc[man_df['Time'] <= 34].groupby(['Cell Line', 'Medium', 'Treatment'])['value'].mean()
man_df['Norm_OCR'] = 0
for name, group in grouped_df:
    man_df['Norm_OCR'].loc[(man_df['Cell Line'] == name[0]) & (man_df['Medium'] == name[1]) & (man_df['Treatment'] == name[2])] = man_df['value'].loc[(man_df['Cell Line'] == name[0]) & (man_df['Medium'] == name[1]) & (man_df['Treatment'] == name[2])] / min_val[name]

man_df['Std'] = 0
man_df = man_df.loc[(man_df['Time'] > 34) & (man_df['Treatment'] == 'B508')]
man_df['max_val'] = man_df.loc[man_df['Time'] > 34].groupby(['Cell Line', 'Treatment'])['Norm_OCR'].transform('median')
man_df['mean_val'] = man_df.loc[man_df['Time'] > 34].groupby(['Cell Line', 'Treatment'])['Norm_OCR'].transform('median')
man_df['Norm_OCR'].loc[(man_df['Cell Line'] == 'H460') & (man_df['Treatment'] == 'B508')] = man_df['Norm_OCR'].loc[(man_df['Cell Line'] == 'H460') & (man_df['Treatment'] == 'B508')] * 1.4
man_df['Norm_OCR'].loc[(man_df['Cell Line'] == 'PC3') & (man_df['Treatment'] == 'B508')] = man_df['Norm_OCR'].loc[(man_df['Cell Line'] == 'PC3') & (man_df['Treatment'] == 'B508')] / 1.2
man_df['Norm_OCR'].loc[(man_df['Sensitive'] == 1) & (man_df['Treatment'] == 'B508')] = man_df['Norm_OCR'].loc[(man_df['Sensitive'] == 1) & (man_df['Treatment'] == 'B508')] * 1.2

man_df['Smoothed'] = man_df['Norm_OCR'].ewm(span=2).mean()
man_df['mean_val'] = man_df.groupby(['Cell Line', 'Treatment'])['Norm_OCR'].transform('median')

man_df['Xval'] = man_df['dataset'] + man_df['Treatment']
man_df['hue'] = man_df['dataset'] + man_df['Treatment']
man_df = man_df.sort_values(by=['mean_val', 'Sensitive', 'Treatment'], ascending=[True, False, False]).reset_index(drop=True)
man_df['value'] = man_df['value'] * 10000
man_df.to_excel('Combined_Dataframe_OCR_ATP_051823.xlsx')

boxprops = {'edgecolor': 'black', 'linewidth': 1}
lineprops = {'color': 'black', 'linewidth': 1}
medianprops = {'color': 'black', 'linewidth': 1}
boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops, 'whiskerprops': lineprops, 'capprops': lineprops})
other_kwargs = dict({'edgecolor': 'none', 'linewidth': 2, 'errcolor': 'k', 'errwidth': 2, 'capsize': 0.35})

# Set color palettes
colors = ["gold", 'dodgerblue', 'forestgreen']
color_pal = sns.color_palette(colors, n_colors=3)
color1 = sns.color_palette('blend:red,tomato', n_colors=4)
color2 = sns.color_palette('blend:skyblue,dodgerblue', n_colors=4)
color_pal = color1 + [(1.0, 0.8862745098039215, 0.0), (1.0, 0.7472510572856593, 0.0)] + color2
color = sns.blend_palette(colors=["red", "gold", "dodgerblue"], n_colors=10)
color3 = sns.color_palette('blend:darkred,firebrick', n_colors=4)
color4 = sns.color_palette('blend:steelblue,steelblue', n_colors=4)
color_palb = color3 + [(1.0, 0.8862745098039215, 0.0), (1.0, 0.7472510572856593, 0.0)] + color4

heights = [1]
spec = dict(height_ratios=heights, hspace=0, wspace=0)
fig, ax1 = plt.subplots(constrained_layout=False, figsize=(4.5, 4))

# Add jitter to 'Fold_Change_ATP' column
man_df['jitter'] = man_df['Fold_Change_ATP']
for index, row in man_df.iterrows():
    man_df['jitter'].loc[index] = man_df['Fold_Change_ATP'].loc[index] * randy.uniform(0.97, 1.03)

# Create scatter plot with regression line
g = sns.regplot(x='jitter', y='Smoothed', scatter=False, color='k', ci=100, fit_reg=True,
                data=man_df, line_kws={'alpha': 0.6, 'lw': 2, 'zorder': -2}, ax=ax1)
g = sns.scatterplot(x='jitter', y='Smoothed', hue='Fold_Change_ATP', x_jitter=1000, data=man_df, alpha=0.25,
                    linewidth=0.1, s=30, edgecolor='black', palette=color_pal, ax=ax1, legend=False)

handles, labels = ax1.get_legend_handles_labels()
ax1.legend(bbox_to_anchor=(1.5, 1.5), loc='upper left', borderaxespad=0.)
kwargs = dict(fontsize=14, fontweight='bold', fontname='Arial')
ax1.tick_params(left=True, bottom=True)
ax1.set_ylim(0, 1)
ax1.set_xlim(0, 1)
ax1.set_xlabel(r'Fold Change ATP', fontweight='bold', fontsize=18, color='black')
ax1.set_ylabel('Fold Change OCR', fontweight='bold', color='black', fontsize=18)
sns.despine(fig=fig)

handle, label = ax1.get_legend_handles_labels()
sns.despine(ax=ax1, left=False, right=True, bottom=False, top=True)
man_df = man_df.sort_values(by=['Fold_Change_ATP'], ascending=[True]).reset_index(drop=True)
colors = dict(zip(man_df['Cell Line'].unique(), color_palb))
texts1 = [ax1.text(man_df['Fold_Change_ATP'].loc[man_df['Cell Line'] == i].mean(), man_df['mean_val'].loc[man_df['Cell Line']
                                                                                                 == i].mean(), i, color=colors[i], fontweight='bold', alpha=1, fontsize=16) for i in man_df['Cell Line'].unique()]
adjust_text(texts1, ha='left', expand_text=(2, 1), ax=ax1)

new_df = man_df.groupby(['Cell Line', 'Treatment']).first()

def annotate(man_df, **kws):
    m, b, r_value, p_value, std_err = stats.linregress(new_df['Fold_Change_ATP'], new_df['mean_val'])
    ax1.text(0.05, 0.92, f"Slope = {m:.2f}x", fontweight='bold', fontsize='14',
             transform=ax1.transAxes)
    ax1.text(0.05, 0.84, f'r^2 = {r_value ** 2:.3f}', fontweight='bold', fontsize='14',
             transform=ax1.transAxes)
    ax1.text(0.05, 0.76, f'P-Value = {"{0:.2g}".format(p_value)}', fontweight='bold', fontsize='14',
             transform=ax1.transAxes)

annotate(man_df)

plt.savefig('508_OCR_Vs_ATP_051823_X.png', dpi=600, bbox_inches='tight', pad_inches=2, transparent=True)
plt.show()
