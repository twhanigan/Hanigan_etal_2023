import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from adjustText import adjust_text

# Read the data from Excel
df_final = pd.read_excel('Analyzed_Biochem_072122.xlsx')

# Calculate the average and standard deviation
df_final['Average'] = df_final.groupby(['Replicate', 'Treatment', 'Complex'])['Value'].transform('mean')
df_final['Std_dev'] = df_final.groupby(['Replicate', 'Treatment', 'Complex'])['Value'].transform('mean')

# Define colors for the plot
colors = ["red", 'dodgerblue', 'gold', 'forestgreen', 'orange', 'k']
color_pal = sns.color_palette(colors, n_colors=6)

# Set boxplot properties
boxprops = {'edgecolor': 'black', 'linewidth': 1}
lineprops = {'color': 'black', 'linewidth': 1}
medianprops = {'color': 'black', 'linewidth': 1}
boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops, 'whiskerprops': lineprops, 'capprops': lineprops}
other_kwargs = {'edgecolor': 'k', 'linewidth': 1.5, 'errcolor': 'black', 'errwidth': 2, 'capsize': 0.15,
                'line_kws': {'alpha': 0.1}}

# Define labels and heights for the plot
labels_new = ['508 1 uM', '508 0.3 uM', '508 0.1 uM', '508 0.03 uM', 'Control', 'Known Inhibitor']
heights = [1]

# Create the figure and axes
spec = {'height_ratios': heights, 'hspace': 0.1, 'wspace': 0}
fig, ax1 = plt.subplots(nrows=1, ncols=1, constrained_layout=False, sharex=True, sharey=True, figsize=(4, 4),
                         gridspec_kw=spec)

# Plot the barplot and stripplot
g = sns.barplot(x='Complex', y='Value', hue='Treatment', edgecolor=None, linewidth=1, data=df_final, dodge=True,
                palette=color_pal, ax=ax1, alpha=0.75)
g = sns.stripplot(x='Complex', y='Value', hue='Treatment', edgecolor='k', linewidth=1, data=df_final, dodge=True,
                  palette=color_pal, size=5, jitter=0.2, ax=ax1, alpha=0.9)

# Remove the legend from the plot
ax1.get_legend().remove()
plt.setp([g.lines], alpha=0.5)

# Get the legend handles and labels
handles1, labels = ax1.get_legend_handles_labels()
kwargs = {'fontsize': 14, 'fontweight': 'bold', 'fontname': 'Arial'}

# Set labels and formatting for the axes
ax1.set_xlabel('')
ax1.set_ylabel('Fold Change Activity\n(nmol/min/mg)', fontsize=18, fontname="Arial", fontweight='bold')
ax1.set_xticklabels(['CI', 'CII', 'CIII', 'CIV'], fontsize=18, fontweight='bold', ha='right')

# Function to change the width of the bars
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value
        patch.set_width(new_value)
        patch.set_x(patch.get_x() + diff * 0.5)

# Adjust the width of the bars
change_width(ax1, 0.13)

# Create a separate legend for custom labels
l1 = plt.legend(handles1[6:12], labels[6:12], bbox_to_anchor=(1.0, 0.5), frameon=True, ncol=3)

# Show only the outside spines
all_axes = fig.get_axes()
for ax in all_axes:
    for sp in ax.spines.values():
        sp.set_visible(True)
    if ax.is_first_row():
        ax.spines['bottom'].set_visible(True)
    if ax.is_first_col():
        ax.spines['left'].set_visible(True)

# Save and show the plot
plt.savefig('Analyzed_Biochem_final_072122E.png', dpi=600, transparent=True, bbox_inches='tight', pad_inches=2)
plt.show()
