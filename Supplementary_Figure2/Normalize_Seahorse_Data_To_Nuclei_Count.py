import numpy as np
import pandas as pd
import math
from pandas import ExcelWriter
from pandas import ExcelFile
import warnings
import seaborn as sns
import statsmodels.api as sm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import style

# Set font family and axes linewidth for plots
plt.rcParams['font.family'] = 'Arial'
plt.rc('axes', linewidth=2)

# Read cell count and Seahorse data
cell_count_df = pd.read_csv('Seahorse_Image_082419.csv')
Seahorse_df = pd.read_excel('Seahorse_Data_082419.xlsx')

# Filter out unwanted groups from Seahorse data
Seahorse_df = Seahorse_df.loc[(Seahorse_df['Group'] != 'Group 21') & (
    Seahorse_df['Group'] != 'Background')]

# Group cell count data by well and calculate average count
grouped = cell_count_df.groupby('Metadata_Well')
avg_well = grouped['Count_DNA'].mean()
mean_df = avg_well.reset_index()


def normalize_seahorse(df1, df2):
    cols = df1.columns.values.tolist()
    cols.append('Normalized_OCR')
    cols.append('Normalized_ECAR')
    df1['Normalized_OCR'] = 'NaN'
    df1['Normalized_ECAR'] = 'NaN'
    df1['Cell_Count'] = 'NaN'
    lst = []
    Normalized_Seahorse = pd.DataFrame(columns=cols)
    for well in df1.Well.unique():
        for time in df1.Time.unique():
            val_ocr = float(df1['OCR'].loc[(df1['Well'] == well) & (df1['Time'] == time)].values /
                            df2['Count_DNA'].loc[df2['Metadata_Well'] == well].mean())
            val_ecar = float(df1['ECAR'].loc[(df1['Well'] == well) & (df1['Time'] == time)].values /
                             df2['Count_DNA'].loc[df2['Metadata_Well'] == well].mean())
            df1['Normalized_OCR'].loc[(df1['Well'] == well) & (
                df1['Time'] == time)] = val_ocr
            df1['Normalized_ECAR'].loc[(df1['Well'] == well) & (
                df1['Time'] == time)] = val_ecar
            df1['Cell_Count'].loc[(df1['Well'] == well) & (
                df1['Time'] == time)] = df2['Count_DNA'].loc[df2['Metadata_Well'] == well].mean()
    return df1


# Normalize Seahorse data
Normalized_Seahorse_df = normalize_seahorse(Seahorse_df, mean_df)
print(Normalized_Seahorse_df.head())

# Save normalized Seahorse data to Excel
Normalized_Seahorse_df.to_excel('NormalizedcellNum_082419.xlsx')

# Assign Group B values based on Group column
Normalized_Seahorse_df['Group B'] = 'NaN'
Normalized_Seahorse_df['Group B'].loc[(
    Normalized_Seahorse_df['Group'].str.contains('H460'))] = 1
Normalized_Seahorse_df['Group B'].loc[Normalized_Seahorse_df['Group'].str.contains(
    'PC3')] = 2
Normalized_Seahorse_df['Group B'].loc[Normalized_Seahorse_df['Group'].str.contains(
    'KP4')] = 3
Normalized_Seahorse_df['Group B'].loc[Normalized_Seahorse_df['Group'].str.contains(
    'SW620')] = 4
Normalized_Seahorse_df['Group B'].loc[Normalized_Seahorse_df['Group'].str.contains(
    'H2122')] = 5


def func(x, a, b, c):
    return a * np.exp(-b * x) + c


# Multiply Seahorse data to reflect 1000 cells
Normalized_Seahorse_df['Normalized_OCR'] = Normalized_Seahorse_df['Normalized_OCR']*1000
Normalized_Seahorse_df['Normalized_ECAR'] = Normalized_Seahorse_df['Normalized_ECAR']*1000


# Plot line plots for Normalized_OCR
fig, ax = plt.subplots(constrained_layout=False, figsize=(5, 4))
g = sns.lineplot(x='Time', y='Normalized_OCR', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 1],
                 hue='Group', palette='Reds', markers="o", legend=True, ax=ax)
g = sns.lineplot(x='Time', y='Normalized_OCR', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 2],
                 hue='Group', palette='Blues', markers="o", legend=True, ax=ax)
g = sns.lineplot(x='Time', y='Normalized_OCR', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 3],
                 hue='Group', palette='Reds', markers="o", legend=True, ax=ax)
g = sns.lineplot(x='Time', y='Normalized_OCR', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 4],
                 hue='Group', palette='Blues', markers="o", legend=True, ax=ax)
g = sns.lineplot(x='Time', y='Normalized_OCR', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 5],
                 hue='Group', palette='YlOrBr', markers="o", legend=True, ax=ax)
ax.legend(bbox_to_anchor=(1.5, 1.5), loc='upper left', borderaxespad=0.)
ax.axvline(34, linewidth=0.5, color='black', lw=2, ls='--', alpha=0.5)
ax.text(34, 180, 'Drug Treatment', rotation=45,
        fontsize=18, fontname="Arial", fontweight='bold')
plt.xlabel('Time (min)', fontsize=18, fontname="Arial", fontweight='bold')
ax.set_xticklabels(['','0','10','20','30','40'],fontsize=12)
plt.ylabel('OCR (pmol/min/1K cells)', fontsize=18,
           fontname="Arial", fontweight='bold')
plt.yticks(fontsize=12)
plt.savefig('OCR_NormCellNum_082419.png', dpi=600,
            transparent=True, bbox_inches='tight', pad_inches=2)
plt.show()


# Plot scatter plots for Cell_Count
fig, ax = plt.subplots(constrained_layout=False, figsize=(5, 4))
g = sns.scatterplot(x='Time', y='Cell_Count', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 1],
                    hue='Group', palette='Reds', markers="o", legend=True, ax=ax)
g = sns.scatterplot(x='Time', y='Cell_Count', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 2],
                    hue='Group', palette='Blues', markers="o", legend=True, ax=ax)
g = sns.scatterplot(x='Time', y='Cell_Count', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 3],
                    hue='Group', palette='Reds', markers="o", legend=True, ax=ax)
g = sns.scatterplot(x='Time', y='Cell_Count', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 4],
                    hue='Group', palette='Blues', markers="o", legend=True, ax=ax)
g = sns.scatterplot(x='Time', y='Cell_Count', data=Normalized_Seahorse_df.loc[Normalized_Seahorse_df['Group B'] == 5],
                    hue='Group', palette='YlOrBr', markers="o", legend=True, ax=ax)
plt.ylim(0, 4000)
ax.legend(bbox_to_anchor=(1.5, 1.5), loc='upper left', borderaxespad=0.)
ax.axvline(34, linewidth=0.5, color='black', lw=2,ls='--', alpha=0.5)
ax.text(34, 4100, 'Drug Treatment', rotation=45,
        fontsize=18, fontname="Arial", fontweight='bold')
plt.yticks(fontsize=12)
ax.set_xticklabels(['','0','10','20','30','40'],fontsize=12)
plt.savefig('CellCount_082419.png', dpi=600, transparent=True,
            bbox_inches='tight', pad_inches=2)
plt.xlabel('Time (min)', fontsize=18, fontname="Arial", fontweight='bold')
plt.show()
