import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Read the CSV file
subsetb = pd.read_csv('Appended_Combined_H2AX_Puncta.csv')

# Define color palette
colors = ["forestgreen", "gold", "dodgerblue"]
color_pal = sns.color_palette(colors, n_colors=3)

# Define boxplot properties
boxprops = {'edgecolor': 'black', 'linewidth': 1}
lineprops = {'color': 'black', 'linewidth': 1}
medianprops = {'color': 'black', 'linewidth': 1}
boxplot_kwargs = dict({'boxprops': boxprops, 'medianprops': lineprops,
                       'whiskerprops': lineprops, 'capprops': lineprops})
other_kwargs = dict({'edgecolor': 'w', 'linewidth': 1.5, 'errcolor': 'black',
                     'errwidth': 2, 'capsize': 0.15, 'line_kws': {'alpha': 0.1}})

heights = [1]
spec = dict(height_ratios=heights, hspace=0, wspace=0.1)

# Create subplots
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, constrained_layout=False,
                               sharex=True, sharey=True, figsize=(4, 4.25), gridspec_kw=spec)

# Plot violin and strip plots for the 'resistant' group
g = sns.violinplot(x='Cpd', y='Puncta_Cell', hue='B508',
                   data=subsetb[subsetb['B508'] == 'resistant'],
                   palette='Blues', ax=ax2, inner=None, alpha=0.1)
g = sns.stripplot(x='Cpd', y='Puncta_Cell', hue='dataset',
                  data=subsetb[subsetb['B508'] == 'resistant'],
                  dodge=True, jitter=0.3, size=3, palette='Blues',
                  alpha=0.6, edgecolor='k', linewidth=0.1, ax=ax2)

# Plot violin and strip plots for the 'sensitive' group
g = sns.violinplot(x='Cpd', y='Puncta_Cell', hue='B508',
                   data=subsetb[subsetb['B508'] == 'sensitive'],
                   palette='Reds', ax=ax1, inner=None, alpha=0.1)
g = sns.stripplot(x='Cpd', y='Puncta_Cell', hue='dataset',
                  data=subsetb[subsetb['B508'] == 'sensitive'],
                  dodge=True, jitter=0.3, size=3, palette='Reds',
                  alpha=0.6, edgecolor='k', linewidth=0.1, ax=ax1)

# Remove legends
ax1.get_legend().remove()
ax2.get_legend().remove()

# Set labels and formatting
ax1.set_ylabel('H2AX Puncta/Cell', fontweight='bold', fontsize=18)
ax2.tick_params(left=False)
ax1.tick_params(left=True)
ax1.set_xlabel('', fontweight='bold', fontsize=18, color='red')
ax2.set_xlabel('', fontweight='bold', fontsize=18, color='dodgerblue')
ax2.set_ylabel('')

# Set x-axis tick labels
ax1.set_xticklabels(['Control', 'B508', 'CDDP'], rotation='45', fontsize=14)
ax2.set_xticklabels(['Control', 'B508', 'CDDP'], rotation='45', fontsize=14)

# Set y-axis limits
plt.ylim(-4, 70)

# Create custom legend
handle, label = ax1.get_legend_handles_labels()
handles, labels = ax2.get_legend_handles_labels()
handy = handles[3:4] + handle[3:4]
labels_new = ['PC3', 'H460']
plt.legend(handles=handy, labels=labels_new, ncol=1,
           loc='upper right', bbox_to_anchor=(2.1, 1.1, 1.1, 1.1))

# Remove unnecessary spines
kwargs = dict(fontsize=14, fontweight='bold', fontname='Arial')
all_axes = fig.get_axes()
for ax in all_axes:
    for sp in ax.spines.values():
        sp.set_visible(False)
    if ax.is_first_row():
        ax.spines['bottom'].set_visible(True)
    if ax.is_first_col():
        ax.spines['left'].set_visible(True)

# Save and show the figure
plt.savefig('H2AX_Puncta_Cell_060723_CDDP_A.png', dpi=600, transparent=True,
            bbox_inches='tight', pad_inches=2)
plt.show()
