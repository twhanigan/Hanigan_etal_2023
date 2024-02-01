import pandas as pd
from bioinfokit import visuz
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
from matplotlib import style
import numpy as np
from adjustText import adjust_text
from pylab import *
rc={"fontname":"Arial",'fontsize':24,'fontweight':'bold','lines.linewidth':2}
plt.rcParams.update({"font.family":'sans-serif','font.sans-serif':"Arial",'font.size':12,'axes.linewidth':2,'axes.labelsize':14})#,'axes.labelweight':'bold'})
plt.rc('axes', linewidth=2)
#merged_mutations_res = pd.read_csv('ExperimentalvsControl_DE_mygene.csv',dtype={'PValue':float,'gneomic_pos.chr':str},on_bad_lines='skip')

dfs_cat = pd.read_csv('Resistant_Final_Annotated.csv',dtype={'logFC':float,'PValue':float,'genomic_pos.chr':str})
lister = ['protein_coding','Mt_rRNA','Mt_tRNA','rRNA','processed_pseudogene','unprocessed_pseudogene']
#lister = ['protein_coding']
dfs_cat = dfs_cat[dfs_cat['ensembl.type_of_gene'].isin(lister)]
dfs_cat = dfs_cat.drop_duplicates()
dfs_cat['symbol'].replace({'COX1':'MT-CO1','COX2':'MT-CO2','COX3':'MT-CO3'},inplace=True)
sens = dfs_cat[['symbol', 'PValue', 'logFC', 'genomic_pos.chr',
                'pathway.kegg',  'go.CC', 'summary']].copy()
sens.columns = ['Name', 'pval', 'fold_change', 'chromosome', 'kegg',
                'go.CC', 'summary']
sens['logPValue'] = np.log10(sens['pval'])
sens['absaf'] = abs(sens['fold_change'])
sens = sens.drop_duplicates(subset='pval')
sens = sens.drop_duplicates(subset='fold_change')
sens = sens.drop_duplicates(subset='Name')

print(dfs_cat.loc[dfs_cat['symbol'].str.contains('MTCO3P')])
#subsetsens = sens.loc[(sens['summary'].str.contains(
#    'mitoch', na=False)) | (sens['chromosome'] == 'MT') | (sens['go.BP'].str.contains('mitoch'))]
#subsetsens = sens.loc[(sens['chromosome'] == 'MT') | (sens['go.CC'].str.contains('mitoch'))].append(sens[sens['Name'].str.contains('MTCO3P|MTCO2P')],ignore_index=True)
subsetsens = sens.loc[sens['chromosome'] == 'MT']
subsetsens['Type'] = "sensitive"
subsetsens[['pval', 'fold_change']].astype(float)
subsetsens['Significant'] = (subsetsens['pval'] < 0.01) & ((subsetsens['fold_change'] > 2) | (subsetsens['fold_change'] < -2))

#subsetres_sub = subsetres_sub.append(subsetres[(
#    subsetres['Name'] == 'MT-CO2')], ignore_index=True).drop_duplicates()
subsetsens_sub = subsetsens.loc[subsetsens['Significant'] == True].copy().reset_index()
#subsetsens_sub = subsetsens.copy().reset_index()

#sens['pval'].loc[(sens['Name'] == 'MTIF3')] = sens['pval'].loc[(sens['Name'] == 'MTIF3')]/10
#sens['fold_change'].loc[(sens['Name'] == 'MTIF3')] = sens['fold_change'].loc[(sens['Name'] == 'MTIF3')]*3

#subsetsens_sub = subsetsens_sub.append(sens[(
#    sens['Name'] == 'MTIF3')], ignore_index=True)

#concat = subsetsens
#concat = concat.drop_duplicates(subset='Name')
#concat_sub = concat.loc[concat['Significant'] == True].copy().reset_index(drop=True)
#concat =concat.sort_values(by=['absaf','logPValue'],ascending=[False,True]).reset_index(drop=True)
#concat =concat.sort_values(by=['logPValue','absaf'],ascending=[True,False]).reset_index(drop=True)
#sens =sens.sort_values(by=['logPValue','absaf'],ascending=[True,False]).reset_index(drop=True)
sens =sens.sort_values(by=['absaf'],ascending=[True]).reset_index(drop=True)


# sign_line=True,genenames=({'ENSG00000241429':'EEF1A1',"ENSG00000240036":'MT-CO2',"ENSG00000211459":"MT-RNR1","ENSG00000198938":"MT-CO3","ENSG00000114923":"SLC4A3","ENSG00000214026":"MRPL23","ENSG00000159423":"ALDH4A1"}),show=True)

style.use('seaborn-muted')
colors = ["skyblue", "red" ]
color_pal = sns.color_palette(colors, n_colors=2)
#custom = sns.palplot(sns.xkcd_palette(colors),2)
order = ['sensitive', 'resistant']
# sns.set(style='white')
#sizes = 2**concat['Fold Change\nResistant Vs. Sensitive']

#g = sns.relplot(x='fold_change', y='pval', data=concat,size='sizes', sizes=(1, 200), hue='dataset',
#                alpha=.75, palette=color_pal, legend='brief', aspect=1.2, linewidth=0.1, edgecolor='black', height=5)

heights = [2]
spec = dict(height_ratios=heights,hspace=0.1,wspace=0)
fig,ax1 = plt.subplots(nrows=1, ncols=1,sharex=True, sharey=True,constrained_layout=True,figsize=(3.5,4),gridspec_kw=spec)
#sizes = 2**abs(concat['Diff_Res_Sens'])

g=sns.scatterplot(x='fold_change', y='pval', data=sens,size='fold_change',
                 legend='brief', linewidth=0.5, edgecolor='black',color='gray', sizes=(30,10) ,alpha=.1,ax=ax1)#size_norm=(-5, 5),
g=sns.scatterplot(x='fold_change',y='pval',data= subsetsens_sub,size='fold_change', sizes=(200,50),alpha=.95, color='red',legend='brief',linewidth=0.5,edgecolor='black',ax=ax1)#size_norm=(-5, 5)

plt.yscale('log')
#plt.xscale("linear")
#plt.ylim(1, 0.00000001)
plt.ylim(1, 0.0000000001)
plt.xlim(-5, 5)

#plt.xlim(-2, 2)
#plt.axvline(-0.03, linewidth=1, color='black', ls='--', alpha=0.4)
#plt.axvline(0.03, linewidth=1, color='black', ls='--', alpha=0.4)
handles, labels = plt.gca().get_legend_handles_labels()
#plt.gca().get_legend().remove()

handles[3].set_color('red')
handles_new = handles[3],handles[0] # ,handle3,handle4,handle5,handle6,
#labels_new = ['Resistant','Sensitive','1-Fold','2-Fold','3-Fold','4-Fold','5-Fold']
labels_new = ['Mitochondria', 'Other']
# create a legend only using the items
#axes.get_legend_remove()
plt.legend(handles_new, labels_new, title='',
           bbox_to_anchor=(1.33, 1.12), fontsize=10)

plt.xlabel(r'Log$_{2}$(B508/Ctrl)', fontsize=18, fontweight='bold', fontname="Arial")
plt.xticks([-5,-2.5,0,2.5,5],fontsize=12)

#plt.xticks([-2,-1,0,1,2],fontsize=12)
#x_ticks = np.arange(-1.25, 1.25, 0.5)
#ax1.set_xticklabels([-2,-1,0,1,2])
plt.ylabel('PValue', fontsize=18, fontweight='bold', fontname="Arial")
plt.yticks(fontsize=12)
texts = subsetsens_sub.sort_values(by=['absaf','logPValue'],ascending=[False,True]).reset_index(drop=True).head(5)
texts1 = [plt.text(texts.fold_change[i], texts.pval[i],
                   texts.Name[i], color='black',fontweight='bold', fontsize=14) for i in range(len(texts))]
#texts1 = [plt.text(subsetsens_sub.fold_change[i], subsetsens_sub.pval[i],
#                   subsetsens_sub.Name[i], color='crimson', fontsize=14) for i in range(len(subsetsens_sub))]
#texts2 = [plt.text(subsetsens_sub.fold_change[i], subsetsens_sub.pval[i],
#                   subsetsens_sub.Name[i], color='crimson', fontsize=14) for i in range(len(subsetsens_sub))]

adjust_text(texts1, arrowprops=dict(arrowstyle='->', color='black', alpha=1),ha='left',force_points=2,)#expand_text=(3, 1)
#adjust_text(texts2, arrowprops=dict(arrowstyle='->', color='crimson', alpha=1),expand_text=(1, 1),add_objects=texts1)
all_axes = fig.get_axes()
for ax in all_axes:
    for sp in ax.spines.values():
        #sp.set_visible(True)
        sp.set_visible(False)
    if ax.is_first_row():
        ax.spines['bottom'].set_visible(True)
        ax.tick_params(left=False)
#        #ax.spines['bottom'].set_position(('outward',10))
#        ax.spines['top'].set_visible(True)
#    if ax.is_last_row():
#        ax.spines['bottom'].set_visible(True)
#        #ax.spines['bottom'].set_position(('outward',10))
    if ax.is_first_col():
        ax.spines['left'].set_visible(True)
        #ax.set_ylabel('log PValue',fontdict={'fontname':"Arial",'fontsize':'16','fontweight':'bold'})
        ax.tick_params(left=True)
#        #ax.spines['left'].set_position(('outward',10))
#    if ax.is_last_col():
#       ax.spines['right'].set_visible(True)
#ax1.get_yaxis().set_tick_params(which='minor', size=5)
plt.savefig('Resistant_Test_Pythonpipeline_mtDNAb.png',dpi=600,bbox_inches ='tight',pad_inches=1,transparent=True)
plt.show()
