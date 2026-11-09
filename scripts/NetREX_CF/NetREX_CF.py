import numpy as np
import pandas as pd
import sys
from implicit import als
import math
import cvxopt
from scipy.linalg import fractional_matrix_power, block_diag
import matplotlib
from pylab import *
import argparse
#import mosek
from scipy.sparse import csr_matrix
import copy
import cvxpy as cvx
import anndata
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import os 
import progressbar
import matplotlib.pyplot as plt

import NetREX_CF_class as NRC

os.environ["MKL_NUM_THREADS"] = '1'



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_exp", type=str, required=True)
    parser.add_argument("--tf_lst", type=str, required=True)
    parser.add_argument("--prior_path", type=str, required=True)
    parser.add_argument("--out_A", type=str, required=True)
    parser.add_argument("--out_P", type=str, required=True)
    parser.add_argument("--out_confidence", type=str, required=True)
    parser.add_argument("--weights", type=int, default=60)

    args = parser.parse_args()

    in_exp=args.in_exp#"/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/inputtoTFA/0.2/k20/Calb/rep0/E.txt"
    tf_lst=args.tf_lst#"/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt"
    prior_path=args.prior_path#"/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/inputtoTFA/0.2/k20/Calb_prior.txt"
    weights=args.weights # as suggested

    raw_dat = pd.read_csv(in_exp, sep='\t', header=None, index_col=0)
    TF_used = np.loadtxt(tf_lst, dtype=str)
    prior_tab = np.loadtxt(prior_path, dtype=str)

    gene_dict={raw_dat.index[i]:i for i in range(raw_dat.shape[0])}
    reg_dict={TF_used[i]:i for i in range(len(TF_used))}

    prior_A = np.zeros((len(gene_dict),len(reg_dict)))
    for r, t, _ in prior_tab:
            if r in reg_dict.keys() and t in gene_dict.keys():
                prior_A[gene_dict[t], reg_dict[r]]=_

    prior_A *= weights

    C_Input = prior_A

    E_Input = raw_dat

    P = (C_Input != 0).astype(float) 

    # Initi Lambda
    Lambda = np.ones(C_Input.shape)*5

    a1 = 10*C_Input
    a2 = np.max(C_Input)#[100]
    a3 = C_Input
    a4 = 1

    reg = 1

    alpha = a3*P+a4*(1-P)
    beta = (a1-a3)*P+(a2-a4)*(1-P)

    Model = NRC.CFNetREX(E_Input.to_numpy(), C_Input, Lambda, alpha, beta, reg, lambda_a = 1., lambda_s = 1., samplesize=100, eta=0.0000001, iternum=100, p= 150)
    Model.NetREX_Plus_algorithm()

    RMat = Model.X.T.dot(Model.Y)

    np.savetxt(args.out_A, Model.S, fmt='%1.5f', delimiter='\t')
    np.savetxt(args.out_P, Model.A, fmt='%1.5f', delimiter='\t')
    np.savetxt(args.out_confidence, RMat, fmt='%1.5f', delimiter='\t')
