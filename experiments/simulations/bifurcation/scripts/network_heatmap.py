import numpy as np
import seaborn as sns
import argparse
import pandas as pd
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
from scipy.stats import entropy

parser = argparse.ArgumentParser()

parser.add_argument("--ingraph", help="input network, a two column file where each line is regulator and target", type=str, required=True)
parser.add_argument("--outdir", help="output directory containing the heatmaps for transcription factors (raw tf by gene graph, nmf, and pca).", type=str, required=True)
parser.add_argument("--prefix", help="prefix for heatmap title", type=str, required=True)

args = parser.parse_args()
prefix = args.prefix
outdir = args.outdir

ingraph = np.loadtxt(args.ingraph, dtype=str)

tf = np.unique(ingraph[:,0])
tgt = np.unique(ingraph[:,1])

id_to_tf = {i:tf[i] for i in range(len(tf))}
tf_to_id = {tf[i]:i for i in range(len(tf))}
id_to_tgt = {i:tgt[i] for i in range(len(tgt))}
tgt_to_id = {tgt[i]:i for i in range(len(tgt))}

df = np.zeros((len(tf), len(tgt)))

for reg, target in ingraph[:,:2]:
    df[tf_to_id[reg],tgt_to_id[target]] = 1

pd_df = pd.DataFrame(df, index=tf, columns=tgt)

out_hm=sns.clustermap(pd_df,xticklabels=False)
out_hm.ax_col_dendrogram.set_visible(False)
out_hm.fig.suptitle(prefix+' tf by gene graph', fontsize=20) 

out_hm.figure.savefig(outdir+'/'+prefix+'_raw.png')

plt.clf()

# inputtoTFA/0.2/k10/Sbay_prior.txt
# ../062024_lda_formulation/data/1000_n_300_p_3000_m_20_k1_10_k2_symm_dirich_5_alpha_2_pi/priors.txt
# testoutput/0.2/Sbay/nw.tsv

##### NMF ######
model = NMF(n_components=10, init='random', random_state=0)
W = model.fit_transform(pd_df)
H = model.components_

W_df = pd.DataFrame(W, index=tf, columns=[str(_) for _ in range(W.shape[1])])

out_hm=sns.clustermap(W_df,xticklabels=False)
out_hm.ax_col_dendrogram.set_visible(False)
out_hm.fig.suptitle(prefix+' nmf', fontsize=20) 

out_hm.figure.savefig(outdir+'/'+prefix+'_nmf.png')
plt.clf()

U, S, Vh = np.linalg.svd(pd_df, full_matrices=True)
U_pd = pd.DataFrame(U, index=tf, columns=[str(_) for _ in range(U.shape[1])])

out_pca=sns.clustermap(U_pd.iloc[:,:5],xticklabels=False)
out_pca.ax_col_dendrogram.set_visible(False)
out_pca.fig.suptitle(prefix+' pca', fontsize=20) 
plt.savefig(outdir+'/'+prefix+'_pca.png')
plt.clf()