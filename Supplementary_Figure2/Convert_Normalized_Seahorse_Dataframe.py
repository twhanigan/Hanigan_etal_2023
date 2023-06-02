import numpy as np
import pandas as pd
import math
from pandas import ExcelWriter
from pandas import ExcelFile
import warnings
import numbers

# Set font family and axes linewidth for plots
plt.rcParams['font.family'] = 'Arial'
plt.rc('axes', linewidth=2)

# Read OCR and ECAR data from Excel files
OCR_df = pd.read_excel('NormalizedcellNum_082519_c.xlsx').drop('Unnamed: 0', axis=1)
ECAR_df = pd.read_excel('NormalizedcellNum_082419_c.xlsx').drop('Unnamed: 0', axis=1)

# Select relevant columns from OCR and ECAR dataframes
OCR = OCR_df[['Group', 'Time', 'OCR', 'ECAR', 'Normalized_OCR', 'Normalized_ECAR']]
ECAR = ECAR_df[['Group', 'Time', 'OCR', 'ECAR', 'Normalized_OCR', 'Normalized_ECAR']]

# Concatenate OCR and ECAR dataframes
concat = pd.concat([OCR, ECAR])
concat = concat.set_index(['Group', 'Time'])

# Define sensitive and resistant cell lines
sensitive = ['H460', 'A549', 'H1975', 'HCT116', 'H2122']
resistant = ['SW620', 'H2009', 'MDAMB231', 'PC3', 'KP4']

# Stack the concatenated dataframe
stacked = concat.stack().reset_index()
stacked.columns = ['Cell Line', 'Time', 'dataset', 'value']

# Extract treatment information from Cell Line column
stacked['Treatment'] = stacked['Cell Line'].str.split('\s', expand=True, n=1)[1]
stacked['Medium'] = 'Basal'
stacked['Medium'].loc[stacked.Treatment.str.contains('gln')] = 'Gln'

# Extract replicate information from Treatment column
stacked['Replicate'] = stacked['Treatment'].str.split('_', expand=True, n=1)[1]
stacked['Treatment'].loc[stacked.Treatment.str.contains('B508')] = 'B508'
stacked['Treatment'].loc[stacked.Treatment.str.contains('Control')] = 'Control'

# Extract cell line information from Cell Line column
stacked['Cell Line'] = stacked['Cell Line'].str.split('\s').str[0]

# Set the Sensitive column based on the cell line
stacked['Sensitive'] = 0
stacked['Sensitive'].loc[stacked['Cell Line'].isin(sensitive)] = 1
stacked['Sensitive'].loc[stacked['Cell Line'].isin(resistant)] = 0

# Remove rows with Cell Line 'Group'
stacked = stacked.loc[stacked['Cell Line'] != 'Group']

# Sort the stacked dataframe
stacked = stacked.sort_values(by=['Cell Line', 'dataset', 'Medium', 'Treatment'])

# Save the stacked dataframe to an Excel file
stacked.to_excel('Combined_RawData_060223_percell_a.xlsx')
