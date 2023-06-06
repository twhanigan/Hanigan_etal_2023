import math
import numpy as np
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import warnings
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Patch
import matplotlib_tools as mpt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox, TextArea
import re
import decimal
from adjustText import adjust_text
from decimal import Decimal
from pylab import *
import rdkit as rd
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdRGroupDecomposition
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import rdMolDraw2D, MolDrawOptions
from cairosvg import svg2png
from PIL import Image

# Define matplotlib plotting properties

rc = {"fontname": "Arial", 'fontsize': 24,
      'fontweight': 'bold', 'lines.linewidth': 2}
plt.rcParams.update({"font.family": 'sans-serif', 'font.sans-serif': "Arial", 'font.size': 12,
                     'axes.linewidth': 2, 'axes.labelsize': 14, 'axes.labelweight': 'bold', 'legend.handlelength': 1})
plt.rc('axes', linewidth=2)
sns.set_style(style="white")

# Import the data table containing the SAR information

df = pd.read_excel('SAR_042823.xlsx')

# Add an Rdkit molecule column

PandasTools.AddMoleculeColumnToFrame(
    df, smilesCol="SMILES", includeFingerprints=True)

# Convert pGC50 to GC50, subset the dataframe and define active compounds in 'color' column

df['GC50'] = np.log10(df['GC50 (H460)'])
df = df.fillna(0)
df_subset = df[['Name', 'GC50', 'Structure.1',
                'GC50 (H460)', 'ROMol', 'SMILES']].copy()
df_new = df_subset.set_index(['Name'])
df_subset['color'] = 'Inactive'
df_subset['color'].loc[(df_subset['GC50'] < -6)] = 'Active'

# convert molecule column to png image

def mol_to_img(mol):
    print(mol.GetBonds())
    AllChem.Compute2DCoords(mol)
    #mol_id = mol.GetProp('SMILES')
    d = rdMolDraw2D.MolDraw2DCairo(1, 1)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    img = d.GetDrawingText()
    return img

# Function to align the molecules to the

def align_molecules(mol, temp):
    AllChem.Compute2DCoords(temp)
    AllChem.Compute2DCoords(mol)
    Aligned = AllChem.GenerateDepictionMatching2DStructure(mol, temp)
    img = Draw.MolToImage(Aligned, dpi=1200)
    return img

# function to convert molecule column to svg image

def moltosvg(mol, molSize=(300, 300), kekulize=True):
    #mc = Chem.Mol(mol.ToBinary())
    mc = mol
    AllChem.EmbedMolecule(mol)
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    img = drawer.GetDrawingText()
    return img

# function to align and adjust bond length of molecules to a template

def adjust_bondlength(mol, temp):
    AllChem.Compute2DCoords(mol)
    AllChem.Compute2DCoords(temp)
    ref_bond = mol.GetBonds()[0]
    ref_length = rdMolTransforms.GetBondLength(
        mol.GetConformer(), ref_bond.GetBeginAtomIdx(), ref_bond.GetEndAtomIdx())
    prob_bond = temp.GetBonds()[0]
    prob_length = rdMolTransforms.GetBondLength(
        temp.GetConformer(), prob_bond.GetBeginAtomIdx(), prob_bond.GetEndAtomIdx())
    ratio = ref_length / prob_length
    print(ref_length, ratio)
    AllChem.GenerateDepictionMatching2DStructure(mol, temp)
    return mol

# function to convert molecule images to RGBA based colors and adjust the background transparency

def white_to_transparency(img):
    img = img.convert("RGBA")
    datas = img.getdata()
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    img.putdata(newData)
    return(img)

# Extract molecules from smiles column of dataframe

mols = [Chem.MolFromSmiles(
    i) for i in df_subset['SMILES'].loc[df_subset['color'] == 'Active']]

# define the template molecule to align structures to

template = Chem.MolFromSmiles('C12=CC=CC=C1CNCC2')

# align and adjust the bondlength of all molecules to the template

molecules = [adjust_bondlength(i, template) for i in mols]


# convert molecules to .png image and exctact as a list

image_list = [Draw.MolToImage(i) for i in molecules]

# extract the GC50 values from active compounds

bartop_labels = df_subset['GC50'].loc[df_subset['color'] == 'Active']

# make the figure and plot the data

colors = ['red', 'k']
color_pal = sns.color_palette(colors, n_colors=3)
fig, ax = plt.subplots(1, 1, sharex=True, figsize=(5, 2.5))
fig.subplots_adjust(hspace=-.7)  # adjust space between axes
g = sns.barplot(x=df.Name, y='GC50', edgecolor='k', linewidth=0.5, data=df_subset, dodge=False,
                ax=ax, hue='color', palette=color_pal, alpha=0.85)  # **kwargs)

# Invert and set the y axis limt

plt.gca().invert_yaxis()
ax.set_ylim(-4, -9)

# annotate the bars with images and adjust the postion of the images on the plot

xticks_locations = []
bartop_label_pattern = '',
xtick_labels = [],
orientation = 'v',
img_scale = 0.125
x_offset = 0
y_offset = 0
bartop_label_rotation = 0
x_tail = 0.1
y_tail = 0.5
x_head = 0.9
y_head = 0.8
dx = x_head - x_tail
dy = y_head - y_tail
frameon = True
n = 0
kwargs = dict({'antiliased': True})
for bar, img, top_label in zip(ax.patches, image_list, bartop_labels):
    n = (n+1)
    bar_x_pos = bar.get_x()
    bar_y_pos = bar.get_y()
    _x_offset = bar.get_x()+bar.get_x()*30
    _y_offset = -1
    _x_act = -7.5+(n/7)
    _x = bar_x_pos + _x_offset
    _y = -30+_y_offset*(_x_offset/1.85)
    xticks_locations.append(_x)
    image = white_to_transparency(img)
    imagebox = OffsetImage(image, zoom=0.25, zorder=10)
    imagebox.image.axes = ax
    ab = AnnotationBbox(imagebox, (bar_x_pos, _x_act), xybox=(_x_offset, _y), xycoords='data', boxcoords='offset points',
                        frameon=False, arrowprops=dict(color='k', arrowstyle='->,head_width=0.035,head_length=0.035', zorder=1, lw=1),
                        box_alignment=(0, -1))
    ax.add_artist(ab)


# hide the spines between ax and ax2

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# don't put tick labels at the top

ax.xaxis.tick_top()
ax.tick_params(labeltop=False)
ax.xaxis.tick_bottom()

# make the legend

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles, labels, title=False, bbox_to_anchor=(1.5, 1.5))

# label and format the x and y axes

plt.ylabel(r'Log GC$_{50}$', fontsize=18, fontname="Arial", fontweight='bold')
plt.xlabel('Compound', fontsize=18, fontname="Arial",
           fontweight='bold', labelpad=20)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

# Save the figure

fig.savefig('SAR_waterfall_051723.png', dpi=600,
            transparent=True, bbox_inches='tight', pad_inches=2)
plt.show()
