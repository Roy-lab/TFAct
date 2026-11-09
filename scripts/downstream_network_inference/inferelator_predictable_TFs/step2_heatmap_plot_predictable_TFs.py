# pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import os
import itertools
import numpy as np
import pandas as pd
from scipy import stats

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-t", "--top", help="top edges picked", type=int, required=True)

args = parser.parse_args()
top = str(args.top)

predTFs = pd.read_csv("inferelator_predictable_tf_nums_top" + top + ".txt", sep='\t', header=None)
predTFs.columns = ['algorithm', 'dataset', 'cell', 'GS_src', 'n_TFs']

ALG_FIX = {
    'ibpg': 'IBPG-p',
    'ibpg_2fct': 'TFAct-F_500_100',
    'ibpg_2fct_100_1000': 'TFAct-F_100_1000',
    'ibpg_2fct_2000_100': 'TFAct-F_2000_100',
    'ibpg_2fct_2000_1000': 'TFAct-F_2000_1000',
    'kl_nmf_2fct':'TFAct-KL_500_100', 
    'kl_nmf_2fct_2000_1000': 'TFAct-KL_2000_1000',
    'regularized_nca_tfa_and_exp': 'regularized_nca'
}

all_outdir="GRN_outputs/algos/analysis/inferelator_predictable_TFs/top_" + top
os.system('mkdir -p ' + all_outdir)

scores=predTFs.copy()
scores.algorithm=scores.algorithm.replace(ALG_FIX)

scores.GS_src = scores.GS_src.replace({'Cus_KO_0.01_Cus_ChIP_2_intersect.gs_network': 'Cus_KO_0.01_Cus_ChIP_2_intersect', 
                                       'Cus_ChIP_2.gs_network': 'Cus_ChIP_2',
                                       'Cus_KO_0.01.gs_network': 'Cus_KO_0.01',
                                       'Cus_KO_0.01_union_intersect.gs_network':'Cus_KO_0.01_union_intersect',
                                       'UniBind.gs_network':'UniBind',
                                       "mESC_chipunion_KDUnion_intersect":"ChIP_KD_Intersect",
                                       'mESC_KDUnion':'KD',
                                       'mESC_chipunion':'ChIP'})


scores = scores.groupby(['algorithm', 'dataset', 'GS_src'], as_index=False)['n_TFs'].sum()

METRICS = [
    'n_TFs'
]
algorithm_orders=['default', 'no_TFA', 'regularized_nca', 'regularized_nca_tfa_only', 'TFAct-F_500_100', 'TFAct-F_2000_1000', 'TFAct-F_2000_100', 'TFAct-F_100_1000', 'TFAct-KL_500_100', 'TFAct-KL_2000_1000', 'NINCA_quad', 'TIGER', 'supirfactor_hierarchical', 'supirfactor_shallow', 'NetREX_CF']

ALG_PALETTE = sns.color_palette('colorblind', n_colors=len(algorithm_orders))
PALETTE_MAP = {a: c for a, c in zip(algorithm_orders, ALG_PALETTE)}

scores = pd.melt(scores, 
                 id_vars=['dataset', 'algorithm', 'GS_src'], # The columns conserved
                 value_vars=METRICS, # The columns to melt
                 var_name='metric', # The summary column name of the category
                 value_name='score')


scores['rank'] = scores.groupby(['dataset', 'GS_src', 'metric'])['score']\
                       .rank(ascending=False, method='average')
human_PRIMARY_GS = ['Cus_ChIP_2', 'Cus_KO_0.01',
       'Cus_KO_0.01_Cus_ChIP_2_intersect', 'Cus_KO_0.01_union_intersect', 'UniBind']
mouse_PRIMARY_GS = ['KD', 'ChIP', 'ChIP_KD_Intersect']
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

OUTDIR=all_outdir
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

OUTDIR=all_outdir
plt.savefig(os.path.join(OUTDIR, 'mouse_num.pdf'), bbox_inches='tight')



###
## Human
###
human_DATASETS=['BH', 'RFH']
human_ROW_INDEX = ['Cus_ChIP_2', 'Cus_KO_0.01', 'Cus_KO_0.01_Cus_ChIP_2_intersect', 'Cus_KO_0.01_union_intersect', 'UniBind']

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
mouse_ROW_INDEX = ['ChIP', 'KD', 'ChIP_KD_Intersect']

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