import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import cycle

# Set matplotlib parameters
rc = {"fontname": "Arial", "fontsize": 24, "fontweight": "bold", "lines.linewidth": 2}
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": "Arial", "font.size": 12, "axes.linewidth": 2, "axes.labelsize": 14})
plt.rc("axes", linewidth=2)

Total_Final = pd.read_csv("CMTMRos_Integrated_OXA1LVar.csv")
#list1 = ["PC3", "Parental H460"]
list1 = ['OXA1L S423S','OXA1L S423InsAGC',]
Total_Final = Total_Final.loc[Total_Final["hue"].isin(list1)]
row_nums = ["Ctrl", "10 uM B508", "3 uM B508", "1 uM B508", "0.3 uM B508", "0.1 uM B508"]
row_cycle = cycle(row_nums)
print(row_nums)
print(next(row_cycle))

Concentrations = ["1 uM B508", "0.3 uM B508", "0.1 uM B508", "0.03 uM B508", "0.01 uM B508", "Ctrl"]
for_order = ["3 \u03BCM", "1 \u03BCM", "0.3 \u03BCM", "0.1 \u03BCM", "0.03 \u03BCM", "Ctrl"]
rev_order = ["Ctrl", "0.1 \u03BCM", "0.3 \u03BCM", "1 \u03BCM", "3 \u03BCM", "10 \u03BCM"]

Conc = [0.00001, 3.3333e-6, 1.11111e-6, 3.7037e-7, 1.2345e-7, 4.1152e-9]
LogConc = np.log10(Conc)

colors = ["tomato", "dodgerblue"]
colorsb = ["red", "dodgerblue"]

color_pal = sns.color_palette(colors, n_colors=2)
color_pal2 = sns.color_palette(colorsb, n_colors=2)


def dist_plot_grid(frame, name):
    g = sns.FacetGrid(
        frame,
        row="Concentration",
        col="hue",
        hue="dataframe",
        col_order=list1,
        sharex=True,
        sharey=True,
        aspect=2.45,
        height=0.65,
        palette=color_pal,
        margin_titles=True,
    )
    g.map(
        sns.kdeplot,
        "Intensity",
        label=False,
        bw_adjust=0.2,
        fill=True,
        linewidth=1,
        alpha=0.3,
        log_scale=True,
    )
    g.map(
        sns.kdeplot,
        "Intensity",
        label=False,
        bw_adjust=0.2,
        fill=False,
        linewidth=2,
        alpha=0.7,
        log_scale=True,
    )
    plt.xlim(70, 3000)
    plt.xlabel("Integrated Intensity")
    plt.ylabel("Count")
    g.map(plt.axhline, y=0, lw=2, clip_on=False, color="k", alpha=0.8)
    g.despine(bottom=True, left=True)
    g.set(xticks=[10 ** 2, 10 ** 3], yticks=[1, 2], yticklabels=["", "", ""], ylabel=None, xlabel=None)
    num = -1
    alt = -1
    for ax in g.axes.flat:
        if ax.is_last_col():
            num = num + 1
            ax.text(
                0.9,
                0.25,
                for_order[num],
                fontweight="bold",
                ha="left",
                va="center",
                transform=ax.transAxes,
            )
        if ax.is_first_col():
            ax.spines["left"].set_visible(True)
    for ax in g.axes.flat:
        alt += 1
        if alt % 2 == 0:
            print(alt)
            f = next(row_cycle)
            ax.axvline(
                x=Total_Final["Intensity"].loc[
                    (Total_Final["Concentration"] == f)
                    & (Total_Final["dataframe"] == "B508")
                    & (Total_Final["hue"] == "PC3")
                ].median(),
                color="blue",
                lw=3,
                ls="solid",
                alpha=0.95,
                zorder=-5,
            )
            ax.axvline(
                x=Total_Final["Intensity"].loc[
                    (Total_Final["Concentration"] == "Ctrl") & (Total_Final["hue"] == "PC3")
                ].median(),
                color="red",
                lw=3,
                ls="solid",
                alpha=0.95,
                zorder=-5,
            )
        else:
            ax.axvline(
                x=Total_Final["Intensity"].loc[
                    (Total_Final["Concentration"] == f)
                    & (Total_Final["dataframe"] == "B508")
                    & (Total_Final["hue"] == "Parental H460")
                ].median(),
                color="blue",
                lw=3,
                ls="solid",
                alpha=0.95,
                zorder=-5,
            )
            ax.axvline(
                x=Total_Final["Intensity"].loc[
                    (Total_Final["Concentration"] == "Ctrl") & (Total_Final["hue"] == "Parental H460")
                ].median(),
                color="red",
                lw=3,
                ls="solid",
                alpha=0.95,
                zorder=-5,
            )
    g.fig.subplots_adjust(hspace=-0.4, wspace=-0.2)
    g.set(facecolor="w", alpha=0.1)
    g.set_titles(col_template="", row_template="", size=16)
	g.axes[0, 0].set_title(r'S419(AGC)$_{4}$    ', color='steelblue', fontdict={
	                       'fontsize': 16, 'fontweight': 'bold', 'fontname': "Arial"}, pad=20)  # "rotation":45
	g.axes[0, 1].set_title(r'   S419(AGC)$_{5}$', color='red', fontdict={
	                       'fontsize': 16, 'fontweight': 'bold', 'fontname': "Arial"}, pad=20)
    g.add_legend(title="Compound")
    g.axes[3, 0].set_ylabel("   Cell Count", labelpad=10, fontsize=18, fontname="Arial", fontweight="bold")
    g.fig.text(
        0.42, -0.07, "Log CMTMRos Intensity", ha="center", va="center", fontsize=18, fontname="Arial", fontweight="bold"
    )
    plt.savefig("CMTMRos_" + name + ".png", dpi=600, transparent=True, bbox_inches="tight", pad_inches=1)


dist_plot_grid(Total_Final, "OXA1L_Variants")
