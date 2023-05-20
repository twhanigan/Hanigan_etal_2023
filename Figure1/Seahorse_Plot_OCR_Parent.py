import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read Excel data into a pandas DataFrame
Seahorse_df = pd.read_excel('TWH_B819_OXA1L4S_OXA1L5S_Parent_Dose_110421B.xlsx', sheet_name='Rate')

# Filter out unwanted groups from the DataFrame
Seahorse_df = Seahorse_df.loc[~Seahorse_df['Group'].isin(['Group 21', 'Background', 'f'])]

# Assign numerical labels to groups in 'Group B' column
Seahorse_df['Group B'] = 'NaN'
Seahorse_df.loc[Seahorse_df['Group'].str.contains('H460 OXA1L 4S'), 'Group B'] = 1
Seahorse_df.loc[Seahorse_df['Group'].str.contains('H460 OXA1L 5S'), 'Group B'] = 2
Seahorse_df.loc[Seahorse_df['Group'].str.contains('H460 Parental'), 'Group B'] = 3

# Convert 'Group' column to string type
Seahorse_df['Group'] = Seahorse_df['Group'].astype(str)

# Group the DataFrame by 'Group'
grouped_df = Seahorse_df.groupby('Group')

# Find the minimum and maximum values of OCR for each group
min_val = Seahorse_df.loc[Seahorse_df['Time'] == 79.9205907233333].groupby('Group').OCR.mean()
max_val = Seahorse_df.loc[Seahorse_df['Time'] == 1.30624229333333].groupby('Group').OCR.mean()
print(min_val, max_val)

# Calculate normalized values and smoothed values
Seahorse_df['Norm'] = Seahorse_df.groupby('Group')['OCR'].transform(lambda x: x / max_val[x.name])
Seahorse_df['Smoothed'] = Seahorse_df.groupby('Group')['Norm'].transform(lambda x: x.ewm(span=2).mean())

# Calculate standard deviation, mean, upper and lower bounds for each time and group
Seahorse_df['Std'] = Seahorse_df.groupby(['Time', 'Group'])['OCR'].transform('std')
Seahorse_df['Mean'] = Seahorse_df.groupby(['Time', 'Group'])['Norm'].transform('mean')
Seahorse_df['Upper'] = Seahorse_df['Mean'] + 0.1
Seahorse_df['Lower'] = Seahorse_df['Mean'] - 0.1

# Create the plot
fig, ax = plt.subplots(figsize=(5, 4))
hue_order = list(Seahorse_df['Group'].loc[Seahorse_df['Group B'] == 2].unique())[::-1]
hue_orderb = list(Seahorse_df['Group'].loc[Seahorse_df['Group B'] == 1].unique())[::-1]

# Plot the smoothed values for Group B = 1 and Group B = 2
sns.lineplot(x='Time', y='Smoothed', data=Seahorse_df.loc[Seahorse_df['Group B'] == 1],
             hue='Group', palette='Blues', markers=True, ci=None, marker='o', mec='k',
             estimator=np.mean, err_kws={'alpha': 0.03}, hue_order=hue_orderb, lw=3, ax=ax, markersize=5)

sns.lineplot(x='Time', y='Smoothed', data=Seahorse_df.loc[Seahorse_df['Group B'] == 2],
             hue='Group', palette='Reds', markers=True, ci=None, marker='o', mec='k', estimator=np.mean,
             err_kws={'alpha': 0.03}, hue_order=hue_order, lw=3, ax=ax, markersize=5)

# Creating color palettes
colors = dict(zip(Seahorse_df['Group'].unique(), 'Reds'))
length = len(Seahorse_df['Group'].loc[Seahorse_df['Group B'] == 2].unique())
lengthb = len(Seahorse_df['Group'].loc[Seahorse_df['Group B'] == 1].unique())
m_cmap = sns.color_palette("Reds", length)
m_cmapb = sns.color_palette("Blues", lengthb)

n = -1
groups = list(Seahorse_df['Group'].loc[(Seahorse_df['Group B'] == 2)].unique())
for name, group in Seahorse_df.loc[Seahorse_df['Group B'] == 2].groupby('Group'):
    n = n + 1
    plt.fill_between(Seahorse_df['Time'].loc[Seahorse_df['Group'] == name], Seahorse_df['Lower'].loc[Seahorse_df['Group'] == name], Seahorse_df['Upper'].loc[Seahorse_df['Group'] == name], color=m_cmap[n], alpha=.1)

n = -1
for name, group in Seahorse_df.loc[Seahorse_df['Group B'] == 1].groupby('Group'):
    n = n + 1
    plt.fill_between(Seahorse_df['Time'].loc[Seahorse_df['Group'] == name], Seahorse_df['Lower'].loc[Seahorse_df['Group'] == name], Seahorse_df['Upper'].loc[Seahorse_df['Group'] == name], color=m_cmapb[n], alpha=.1)

# Remove legend
ax.get_legend().remove()

# Vertical lines
plt.axvline(x=14, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)
plt.axvline(x=30, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)
plt.axvline(x=45, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)
plt.axvline(x=63, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)

# Text annotations
plt.text(14.1, 2.4, 'B508', rotation=45, fontsize=18, fontname="Arial", fontweight='bold')
plt.text(30.1, 2.4, 'Oligomycin', rotation=45, fontsize=18, fontname="Arial", fontweight='bold')
plt.text(45.1, 2.4, 'FCCP', rotation=45, fontsize=18, fontname="Arial", fontweight='bold')
plt.text(63.1, 2.4, 'Antimycin\nRotenone', rotation=48, fontsize=16, fontname="Arial", fontweight='bold')

# Legend
handles, labels = ax.get_legend_handles_labels()
labels_new = ['B508 3 µM', 'B508 1 µM', 'B508 0.3 µM', 'B508 0.1 µM', 'B508 0.03 µM', 'B508 0.01 µM', 'B508 0.003 µM', '143-01 3 µM', 'H460 Parental Control']
labels_new = labels_new[::-1]
fig.legend(handles, labels_new, loc='upper right', bbox_to_anchor=(1.1, 1.1, 1.1, 1.1))

# Axis labels and formatting
plt.xlabel('Time (min)', fontsize=18, fontname="Arial", fontweight='bold')
plt.xticks(fontsize=12)
plt.ylabel('Fold Change OCR', fontsize=18, fontname="Arial", fontweight='bold')
plt.yticks(fontsize=12)
plt.ylim(0, 2.3)

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Save and display the plot
plt.savefig('Figure_1f_OCR.png', dpi=600, transparent=True, bbox_inches='tight', pad_inches=2)
plt.show()