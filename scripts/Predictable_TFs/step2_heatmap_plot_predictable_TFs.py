# pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import os
import itertools
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans']

predTFs = pd.read_csv("predictable_tf_nums.txt", sep='\t', header=None)
predTFs.columns = ['algorithm', 'dataset', 'cell', 'GS_src', 'n_TFs']

ALG_FIX = {
    'ibpg': 'IBPG-p',
    'ibpg_2fct': 'TFAct-F_500_100',
    'ibpg_2fct_100_1000': 'TFAct-F_100_1000',
    'ibpg_2fct_2000_100': 'TFAct-F_2000_100',
    'ibpg_2fct_2000_1000': 'TFAct-F_2000_1000',
    'kl_nmf_2fct':'TFAct-KL_500_100', 
    'kl_nmf_2fct_2000_1000': 'TFAct-KL_2000_1000'
}



scores=predTFs.copy()
scores.algorithm=scores.algorithm.replace(ALG_FIX)

scores.GS_src = scores.GS_src.replace({'Cus_KO_0.01_Cus_ChIP_2_intersect.gs_network': 'Cus_KO+Cus_ChIP', 
                                       'Cus_ChIP_2.gs_network': 'Cus_ChIP',
                                       'Cus_KO_0.01.gs_network': 'Cus_KO',
                                       'Cus_KO_0.01_union_intersect.gs_network':'Cus_KO+UniBind',
                                       'UniBind.gs_network':'UniBind',
                                       "mESC_chipunion_KDUnion_intersect":"ChIP+KD",
                                       'mESC_KDUnion':'KD',
                                       'mESC_chipunion':'ChIP'})


scores = scores.groupby(['algorithm', 'dataset', 'GS_src'], as_index=False)['n_TFs'].sum()

METRICS = [
    'n_TFs'
]
algorithm_orders=['HALS-p', 'IBPG-p', 'NINCA_quad', 'NetREX_CF', 'TIGER', 
                  'TFAct-F_100_1000', 'TFAct-F_2000_1000', 'TFAct-F_500_100', 'TFAct-F_2000_100',  
        'TFAct-KL_500_100', 'TFAct-KL_2000_1000',
       'regularized_nca', 'supirfactor_hierarchical', 'supirfactor_shallow']

ALG_PALETTE = sns.color_palette('colorblind') + [sns.xkcd_rgb['dark grey'], sns.xkcd_rgb['teal'],       # Teal
    sns.xkcd_rgb['magenta'],    # Magenta
    sns.xkcd_rgb['gold'] ]
PALETTE_MAP = {a: c for a, c in zip(algorithm_orders, ALG_PALETTE)}

scores = pd.melt(scores, 
                 id_vars=['dataset', 'algorithm', 'GS_src'], # The columns conserved
                 value_vars=METRICS, # The columns to melt
                 var_name='metric', # The summary column name of the category
                 value_name='score')


scores['rank'] = scores.groupby(['dataset', 'GS_src', 'metric'])['score']\
                       .rank(ascending=False, method='average')
human_PRIMARY_GS = ['Cus_ChIP', 'Cus_KO',
       'Cus_KO+Cus_ChIP', 'Cus_KO+UniBind', 'UniBind']
mouse_PRIMARY_GS = ['KD', 'ChIP', 'ChIP+KD']
human_summary = scores.loc[scores.GS_src.isin(human_PRIMARY_GS)]\
                .groupby(['algorithm', 'metric'])['rank'].median().reset_index()

mouse_summary = scores.loc[scores.GS_src.isin(mouse_PRIMARY_GS)]\
                .groupby(['algorithm', 'metric'])['rank'].median().reset_index()
# Re-organize for heatmap
human_summary = human_summary.pivot_table(columns='algorithm', index='metric', values='rank')
human_summary = human_summary[algorithm_orders]

mouse_summary = mouse_summary.pivot_table(columns='algorithm', index='metric', values='rank')
mouse_summary = mouse_summary[algorithm_orders]

human_primary_scores = scores.loc[scores.GS_src.isin(human_PRIMARY_GS)].copy()
mouse_primary_scores = scores.loc[scores.GS_src.isin(mouse_PRIMARY_GS)].copy()


METRIC_TITLES = {
    'n_TFs': 'Predictable TFs',
}

##
# Human
##
fig, ax = plt.subplots(1, 1, figsize=(7, 6))
metric = 'n_TFs'
data = human_primary_scores.loc[(human_primary_scores.algorithm != 'random') &
                              (human_primary_scores.metric == metric)] 
if metric == 'n_TFs':
    yval = 'score'
else:
    yval = 'log_rand_ratio'
    

order = pd.DataFrame(human_summary.iloc[0,:]).median(axis=1).sort_values(ascending=True).index
palette = [PALETTE_MAP[i] for i in order]
        
ax = sns.swarmplot(data=data,
                       x='algorithm', order=order,
                       y=yval,
                       palette=palette,
                       ax=ax)
    
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
[t.set_color(i) for i, t in zip(palette, ax.xaxis.get_ticklabels())]

[t.set_fontweight('bold') for t in ax.xaxis.get_ticklabels()]
#     ax.set_ylabel('Fold change\nrelative to random')
ax.set_ylabel(None)
ax.axhline(0, c='k')
ax.set_xlabel(None)

if metric == 'fscore':
    ax.set_yticks(np.arange(-4, 4, 1))
    ax.set_yticklabels(['1/16', '1/8', '1/4', '1/2', '1', '2', '4', '8'])
    ax.set_ylim(-4, 3)
    ax.set_ylabel('Improvement over random network')
else:
    ax.set_ylabel('# TFs')


ax.set_title(METRIC_TITLES[metric])

OUTDIR="TFA_output/analysis/predictable_TFs/"
plt.savefig(os.path.join(OUTDIR, 'human_num.pdf'), bbox_inches='tight')



##
# mouse
##
fig, ax = plt.subplots(1, 1, figsize=(7, 6))
metric = 'n_TFs'
data = mouse_primary_scores.loc[(mouse_primary_scores.algorithm != 'random') &
                              (mouse_primary_scores.metric == metric)] 
if metric == 'n_TFs':
    yval = 'score'
else:
    yval = 'log_rand_ratio'
    

order = pd.DataFrame(mouse_summary.iloc[0,:]).median(axis=1).sort_values(ascending=True).index
palette = [PALETTE_MAP[i] for i in order]
        
ax = sns.swarmplot(data=data,
                       x='algorithm', order=order,
                       y=yval,
                       palette=palette,
                       ax=ax)
    
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
[t.set_color(i) for i, t in zip(palette, ax.xaxis.get_ticklabels())]

[t.set_fontweight('bold') for t in ax.xaxis.get_ticklabels()]
#     ax.set_ylabel('Fold change\nrelative to random')
ax.set_ylabel(None)
ax.axhline(0, c='k')
ax.set_xlabel(None)

if metric == 'fscore':
    ax.set_yticks(np.arange(-4, 4, 1))
    ax.set_yticklabels(['1/16', '1/8', '1/4', '1/2', '1', '2', '4', '8'])
    ax.set_ylim(-4, 3)
    ax.set_ylabel('Improvement over random network')
else:
    ax.set_ylabel('# TFs')


ax.set_title(METRIC_TITLES[metric])

OUTDIR="TFA_output/analysis/predictable_TFs/"
plt.savefig(os.path.join(OUTDIR, 'mouse_num.pdf'), bbox_inches='tight')



###
## Human
###
human_DATASETS=['BH', 'RFH']
human_ROW_INDEX = ['Cus_ChIP', 'Cus_KO', 'Cus_KO+Cus_ChIP', 'Cus_KO+UniBind', 'UniBind']

fig, axes = plt.subplots(2, 1, figsize=(8, 12))

for i, (dataset, axrow) in enumerate(zip(human_DATASETS, axes)): 
    ax = axrow
    for j, (metric) in enumerate(['n_TFs']):
        df = scores.loc[(scores.dataset == dataset) &
                        (scores.metric == metric) &
                        (scores.GS_src.isin(human_PRIMARY_GS))]
        hmap = df.pivot_table(columns='algorithm', 
                              index='GS_src', 
                              values='rank')
        hmap = hmap.loc[human_ROW_INDEX, order]     
        fmt = '.1f'
        ax = sns.heatmap(hmap, ax=ax,
                         annot=True, fmt=fmt, 
                         cbar=False, mask=hmap.isnull(),
                         cmap='viridis_r', square=True,
                         vmin=0, vmax=11,
                         linewidths=2)
        if j == 0:
            ax.set_yticklabels(human_ROW_INDEX)
            ax.set_ylabel(dataset, rotation=0, ha='right')
        else:
            ax.set_yticklabels([])
            ax.set_ylabel(None)
        if i < 1:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(order, rotation=30, ha='right')
            [t.set_color(PALETTE_MAP[alg]) for alg, t in zip(order, ax.xaxis.get_ticklabels())]
            [t.set_fontweight('bold') for t in ax.xaxis.get_ticklabels()]
        if i == 0:
            ax.set_title(METRIC_TITLES[metric], fontsize=18)
        ax.set_xlabel(None)
        ax.set_facecolor(sns.xkcd_rgb['light grey'])


plt.tight_layout()
plt.subplots_adjust(wspace=0.1)
plt.savefig(os.path.join(OUTDIR, 'humans_rank.pdf'), bbox_inches='tight')


##
## mouse
###
mouse_DATASETS=['FBS']
mouse_ROW_INDEX = ['ChIP', 'KD', 'ChIP+KD']

fig, axes = plt.subplots(2, 1, figsize=(8, 12))

for i, (dataset, axrow) in [(0, ('FBS', axes[0]))]: 
    ax = axrow
    for j, (metric) in enumerate(['n_TFs']):
        df = scores.loc[(scores.dataset == dataset) &
                        (scores.metric == metric) &
                        (scores.GS_src.isin(mouse_PRIMARY_GS))]
        hmap = df.pivot_table(columns='algorithm', 
                              index='GS_src', 
                              values='rank')
        hmap = hmap.loc[mouse_ROW_INDEX, order]     
        fmt = '.1f'
        ax = sns.heatmap(hmap, ax=ax,
                         annot=True, fmt=fmt, 
                         cbar=False, mask=hmap.isnull(),
                         cmap='viridis_r', square=True,
                         vmin=0, vmax=11,
                         linewidths=2)
        if j == 0:
            ax.set_yticklabels(mouse_ROW_INDEX)
            ax.set_ylabel(dataset, rotation=0, ha='right')
        else:
            ax.set_yticklabels([])
            ax.set_ylabel(None)
        if i < 0:
            ax.set_xticklabels([])
        else:
            ax.set_xticklabels(order, rotation=30, ha='right')
            [t.set_color(PALETTE_MAP[alg]) for alg, t in zip(order, ax.xaxis.get_ticklabels())]
            [t.set_fontweight('bold') for t in ax.xaxis.get_ticklabels()]
        if i == 0:
            ax.set_title(METRIC_TITLES[metric], fontsize=18)
        ax.set_xlabel(None)
        ax.set_facecolor(sns.xkcd_rgb['light grey'])


plt.tight_layout()
plt.subplots_adjust(wspace=0.1)
plt.savefig(os.path.join(OUTDIR, 'mouses_rank.pdf'), bbox_inches='tight')