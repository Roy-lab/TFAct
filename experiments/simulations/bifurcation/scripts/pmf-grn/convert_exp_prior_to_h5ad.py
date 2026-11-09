import anndata
import argparse
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--in_exp", help="tab delim expression matrix. feature by sample (cell), without header (cell names).", type=str, required=True)
parser.add_argument("--in_prior", help="tab delim prior edge file with 3 fields: regulator, target, confidence. without header.", type=str, required=True)
parser.add_argument("--cell_names", help="A file where each line is a sample (cell) name. If not specified, cell names will be cell0, cell1, ...", type=str)
parser.add_argument("--tf_lst", help="a column file where each line is the name of a regulator. This specifies the order of tfs in prior network.", type=str, required=True)


parser.add_argument("--out_exp", help="output expression h5ad. sample (cell) by feature. If not supplied then no exp .h5ad will be generated", type=str)
parser.add_argument("--out_prior", help="output prior in h5ad. feature by tf. If not supplied then no prior .h5ad will be generated", type=str)


args = parser.parse_args()

with open(args.in_exp) as f:
    exp = anndata.read_csv(f, delimiter='\t')

if args.cell_names is not None:
    cell_names = np.loadtxt(args.cell_names)
    exp.var_names = cell_names
else:
    exp.var_names = ['cell'+str(_) for _ in range(len(exp.var_names))]


## prior
if args.out_prior is not None:
    tfs = np.loadtxt(args.tf_lst,dtype=str)
    gene = exp.obs_names
    prior_A=np.loadtxt(args.in_prior,dtype=str)
    priorA = np.zeros((exp.shape[0],len(tfs)))

    gene_dict={gene[i]:i for i in range(len(gene))}
    reg_dict={tfs[i]:i for i in range(len(tfs))}

    for r, t, _ in prior_A:
        priorA[gene_dict[t], reg_dict[r]]=_

    prior_ad=anndata.AnnData(priorA)
    prior_ad.obs_names = gene
    prior_ad.var_names = tfs

    prior_ad.write_h5ad(args.out_prior)


# transpose the exp to make it cell by features
exp = exp.transpose()

if args.out_exp is not None:
    exp.write_h5ad(args.out_exp)