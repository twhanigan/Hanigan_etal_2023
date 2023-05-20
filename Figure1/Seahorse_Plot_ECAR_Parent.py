import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read data from Excel file
Seahorse_df = pd.read_excel('TWH_B819_OXA1L4S_OXA1L5S_Parent_Dose_110421B.xlsx', sheet_name='Rate')

# Filter out unwanted groups
Seahorse_df = Seahorse_df.loc[(Seahorse_df['Group'] != 'Group 21') & (Seahorse_df['Group'] != 'Background') & (Seahorse_df['Group'] != 'f')]

# Add 'Group B' column and assign values based on 'Group' column
Seahorse_df['Group B'] = 'NaN'
Seahorse_df['Group B'].loc[Seahorse_df['Group'].str.contains('H460 OXA1L 4S')] = 1
Seahorse_df['Group B'].loc[Seahorse_df['Group'].str.contains('H460 OXA1L 5S')] = 2
Seahorse_df['Group B'].loc[Seahorse_df['Group'].str.contains('H460 Parental')] = 3

# Convert 'Group' column to string type
Seahorse_df.Group = Seahorse_df.Group.astype(str)

# Group data by 'Group'
grouped_df = Seahorse_df.groupby('Group')

# Calculate minimum time value and maximum ECAR value for each group
min_time = Seahorse_df['Time'].min()
max_val = Seahorse_df.loc[Seahorse_df['Time'] == 1.30624229333333].groupby('Group').ECAR.mean()

# Calculate absolute minimum ECAR value for each group
abs_min_val = Seahorse_df.groupby('Group').ECAR.min()

# Add new columns to the DataFrame
Seahorse_df['Norm'] = 0
Seahorse_df['Mean'] = 0
Seahorse_df['Upper'] = 0
Seahorse_df['Lower'] = 0
Seahorse_df['Smoothed'] = 'NaN'

# Calculate normalized, smoothed, standard deviation, mean, upper, and lower values for each group at each time point
for name, group in grouped_df:
    Seahorse_df['Norm'].loc[Seahorse_df['Group'] == name] = (Seahorse_df['ECAR'].loc[Seahorse_df['Group'] == name] / max_val[name])
    Seahorse_df['Smoothed'].loc[Seahorse_df['Group'] == name] = Seahorse_df['Norm'].loc[Seahorse_df['Group'] == name].ewm(span=1).mean()

Seahorse_df['Std'] = 0
for name, group in Seahorse_df.groupby(['Time', 'Group']):
    Seahorse_df['Std'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])] = Seahorse_df['ECAR'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])].std()
    Seahorse_df['Mean'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])] = Seahorse_df['Norm'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])].mean()
    Seahorse_df['Upper'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])] = Seahorse_df['Norm'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])].mean() + 0.1
    Seahorse_df['Lower'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])] = Seahorse_df['Norm'].loc[(Seahorse_df['Group'] == name[1]) & (Seahorse_df['Time'] == name[0])].mean() - 0.1

# Create the plot
fig, ax = plt.subplots(figsize=(5, 4))

# Define hue order for the legend
hue_order = list(Seahorse_df['Group'].loc[(Seahorse_df['Group B'] == 2)].unique())[::-1]
hue_orderb = list(Seahorse_df['Group'].loc[(Seahorse_df['Group B'] == 1)].unique())[::-1]

# Plot lines for each group
sns.lineplot(x='Time', y='Smoothed', data=Seahorse_df.loc[Seahorse_df['Group B'] == 1], hue='Group', palette='Blues', markers=True, ci=None, marker='o', mec='k', estimator=np.mean, err_kws={'alpha': 0.03}, hue_order=hue_orderb, lw=3, ax=ax, markersize=5)
sns.lineplot(x='Time', y='Smoothed', data=Seahorse_df.loc[Seahorse_df['Group B'] == 2], hue='Group', palette='Reds', markers=True, ci=None, marker='o', mec='k', estimator=np.mean, err_kws={'alpha': 0.03}, hue_order=hue_order, lw=3, ax=ax, markersize=5)

# Create color palettes for filling between lines
colors = dict(zip(Seahorse_df['Group'].unique(), 'Reds'))
length = len(Seahorse_df['Group'].loc[Seahorse_df['Group B'] == 2].unique())
lengthb = len(Seahorse_df['Group'].loc[Seahorse_df['Group B'] == 1].unique())
m_cmap = sns.color_palette("Reds", length)
m_cmapb = sns.color_palette("Blues", lengthb)

n = -1
groups = list(Seahorse_df['Group'].loc[(Seahorse_df['Group B'] == 2)].unique())

# Fill between lines for each group
for name, group in Seahorse_df.loc[Seahorse_df['Group B'] == 2].groupby('Group'):
    n = n + 1
    plt.fill_between(Seahorse_df['Time'].loc[Seahorse_df['Group'] == name], Seahorse_df['Lower'].loc[Seahorse_df['Group'] == name], Seahorse_df['Upper'].loc[Seahorse_df['Group'] == name], color=m_cmap[n], alpha=.1)

n = -1
for name, group in Seahorse_df.loc[Seahorse_df['Group B'] == 1].groupby('Group'):
    n = n + 1
    plt.fill_between(Seahorse_df['Time'].loc[Seahorse_df['Group'] == name], Seahorse_df['Lower'].loc[Seahorse_df['Group'] == name], Seahorse_df['Upper'].loc[Seahorse_df['Group'] == name], color=m_cmapb[n], alpha=.1)

ax.get_legend().remove()

# Add vertical lines and text labels
plt.axvline(x=14, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)
plt.axvline(x=30, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)
plt.axvline(x=45, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)
plt.axvline(x=63, color='black', lw=2, ls=':', alpha=0.6, zorder=-2)
plt.text(14.1, 2.4, 'B508', rotation=45, fontsize=18, fontname="Arial", fontweight='bold')
plt.text(30.1, 2.4, 'Oligomycin', rotation=45, fontsize=18, fontname="Arial", fontweight='bold')
plt.text(45.1, 2.4, 'FCCP', rotation=45, fontsize=18, fontname="Arial", fontweight='bold')
plt.text(63.1, 2.4, 'Antimycin\nRotenone', rotation=48, fontsize=16, fontname="Arial", fontweight='bold')

# Adjust the legend labels
handles, labels = ax.get_legend_handles_labels()
print(labels)
labels_new = ['B508 3 μM', 'B508 1 μM', 'B508 0.3 μM', 'B508 0.1 μM', 'B508 0.03 μM', 'B508 0.01 μM', 'B508 0.003 μM', '143-01 3 μM']
labels_new = labels_new[::-2]
fig.legend(handles, labels_new, loc='upper right', bbox_to_anchor=(1.1, 1.1, 1.1, 1.1))

# Set plot labels and limits
plt.xlabel('Time (min)', fontsize=18, fontname="Arial", fontweight='bold')
plt.xticks(fontsize=12)
plt.ylabel('Fold Change ECAR', fontsize=18, fontname="Arial", fontweight='bold')
plt.yticks(fontsize=12)
plt.ylim(0.8, 2.2)

# Display the legend and save the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig('Parent_030723_Final_ECAR_All_H.png', dpi=600, transparent=True, bbox_inches='tight', pad_inches=2)
plt.show()
