import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Import Spreadsheets
Total = pd.read_excel("BrdU_DAPI_Cpd1.xlsx")
Total1 = pd.read_excel("BrdU_DAPI_Cpd2.xlsx")

Concentrations = [
    "10 uM B508",
    "3 uM B508",
    "1 uM B508",
    "0.3 uM B508",
    "0.1 uM B508",
    "0.03 uM B508",
    "0.01 uM B508",
    "0.003 uM B508",
    "0.001 uM B508",
    "0.00033 uM B598",
    "Ctrl",
]
Conc = [
    0.00001,
    3.3333e-6,
    1.11111e-6,
    3.7037e-7,
    1.2345e-7,
    4.1152e-8,
    1.3717e-8,
    4.5724e-9,
    1.5241e-9,
    5.08052e-10,
    1.69e-10,
]
LogConc = np.log10(Conc)
colors = [
    (0.8085351787773933, 0.15501730103806227, 0.1522491349480969),
    (0.9748558246828143, 0.5574009996155325, 0.32272202998846594),
    (0.5564013840830451, 0.7596309111880047, 0.8638985005767013),
    (0.25674740484429065, 0.41522491349480967, 0.6844290657439447),
]

color_pal = sns.color_palette(colors, n_colors=4)


def dist_plot_grid(frame, name):
    # Create a FacetGrid for visualizing subplots
    g = sns.FacetGrid(
        frame,
        row="Concentration",  # Assign 'Concentration' as the row variable
        hue="Concentration",  # Assign 'Concentration' for coloring
        sharex=True,
        sharey=True,
        aspect=1.75,
        height=4,
        palette=color_pal,  # Use the predefined color palette
    )
    
    # Define keyword arguments for scatterplot markers
    kwargs = dict(linewidths=0, edgecolor="none")
    kwargs2 = dict(linewidths=0.1, edgecolor="w")
    
    # Plot the first scatterplot layer using sns.scatterplot()
    g.map(
        sns.scatterplot,
        "Intensity_DAPI",  # X-axis variable
        "Intensity_Brdu",  # Y-axis variable
        palette=color_pal,
        label=False,
        size=2,
        alpha=0.5,
        **kwargs,
    )
    
    # Plot the second scatterplot layer with a white edge using sns.scatterplot()
    g.map(
        sns.scatterplot,
        "Intensity_DAPI",
        "Intensity_Brdu",
        palette=color_pal,
        label=False,
        size=2,
        alpha=0.5,
        **kwargs2,
    )

    # Define a function to label each subplot
    def label(x, color, label):
        ax = plt.gca()  # Get the current axes
        ax.text(
            0.1,
            0.5,
            label,
            fontweight="bold",
            fontname="Arial",
            fontsize=36,
            color="k",
            ha="left",
            va="center",
            transform=ax.transAxes,
        )

    # Apply the labeling function to each subplot
    g.map(label, "Concentration")

    # Adjust the spacing between subplots
    g.fig.subplots_adjust(hspace=-0.25)

    # Set the facecolor and transparency of the grid
    g.set(facecolor="w", alpha=0.1)

    # Remove default titles
    g.set_titles("")

    # Remove y-axis ticks
    g.set(yticks=[])

    # Set x-axis limits
    g.set(xlim=(50, 745))

    # Set y-axis limits
    g.set(ylim=(80, 850))

    # Remove the bottom spine of each subplot
    g.despine(bottom=True)

    # Set x-axis tick positions and labels
    g.set(xticks=[245, 495, 745])
    g.set(xticklabels=["2N", "4N", "6N"])

    # Set x and y-axis labels and formatting
    g.set_axis_labels(
        "Intensity DAPI",  # X-axis label
        "",  # Empty y-axis label
        fontsize=77,
        fontname="Arial",
        labelpad=20,
        fontweight="bold",
    )

    # Add a custom text label to the figure
    g.fig.text(
        0.03,
        0.3,
        "Intensity BrdU",
        fontsize=77,
        fontname="Arial",
        fontweight="bold",
        rotation="vertical",
        color="black",
    )

    # Save the figure
    plt.savefig(
        "Scatter_062722_" + name + "_I.png",
        dpi=700,
        transparent=True,
        bbox_inches="tight",
        pad_inches=2,
    )


dist_plot_grid(Total, "Cpd1")
dist_plot_grid(Total1, "Cpd2")
