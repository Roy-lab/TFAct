import numpy as np
from scipy.stats import truncnorm
import os
import argparse

def int_or_list(value):
    try:
        return int(value)
    except ValueError:
        try:
            return [int(x) for x in value.split(',')]
        except ValueError:
            raise argparse.ArgumentTypeError("Must be an integer or a comma-separated list of integers")

parser = argparse.ArgumentParser()
parser.add_argument("--GRN", help="tab separated 'ground truth', where each line is regulator[tab]target.", type=str, required=True)
parser.add_argument("--gene", help="a column file where each line is the name of a gene. This specifies the order of genes in prior network and output expression.", type=str, required=True)
parser.add_argument("--reg", help="a column file where each line is the name of a regulator. This specifies the order of tfs in prior network.", type=str, required=True)
parser.add_argument("--num_cells", help="number of cells to simulate (default 2000).", type=int, default=2000)
parser.add_argument("--cell_type", help="an integer (default 10; the dirichlet hyperparameters is a vector of 5s of length cell_type) or a list of positive integers separated by comma that serves as dirichlet hyperparameter for cell types. Higher values tend to give less variable cell proportion.", type=int_or_list, default=10)
parser.add_argument("--num_rep", help="an integer specifying the number of TFA replicates (default 10).", type=int, default=10)
parser.add_argument("--outdir", help="output directory that contains ${rep number}/P.txt for tfa and ${rep number}/E.txt for expression.", type=str, required=True)

args = parser.parse_args()

GRN=np.loadtxt(args.GRN, dtype=str, delimiter="\t")
gene=np.loadtxt(args.gene, dtype=str)
reg=np.loadtxt(args.reg, dtype=str)

gene_dict={gene[i]:i for i in range(len(gene))}
reg_dict={reg[i]:i for i in range(len(reg))}

n = len(gene) # matrix A is n by p, where n is the number of genes, and p is the number of transcription factors (regulators)
p = len(reg) # each transcription factor also represents a document of words (genes). (TODO) For now each TF belongs to one topic rather than mixture 
m = args.num_cells # m is the number of cells. The TFA matrix, P, is of dimension (p, m)

if isinstance(args.cell_type, int):
    cell_type=np.array([5] * args.cell_type)
else:
    cell_type=np.array(args.cell_type)

# make the the gene by tf matrix A.
A = np.zeros((n,p))

for r, t in GRN:
    A[gene_dict[t], reg_dict[r]]=1

k2 = len(cell_type) # the number of cell types. We set k2=q, per formulation.jpg within the same directory
epsilon_P = 1 / k2 # the larger the value, the more likely that a tf "belongs" to another cell type get expressed in the current cell type
outdir= args.outdir # change this. Where everything is output

rep = args.num_rep # generate simulated runs

os.system("mkdir -p " + outdir)

pi = [2] * k2 # account for rare cell types. More balanced
rates = [0.1, 0.3, 0.6, 0.9] # drop out expression. Larger values correspond to more dropouts.

for run in range(rep):
    #
    # simulate P
    #
    psi = np.random.dirichlet(pi) 
    c = [np.where(np.random.multinomial(1, psi)==1)[0][0] for i in range(m)] # assign cells to type

    gamma2 = [1/k2] * k2 # TODO [can be balanced or skewed] for assigning tfs to cell types (TODO but each cell has a chance for expressing TF not assigned to its type)
    z_c = np.random.choice(k2, size=p, p=gamma2) # cell type assignment for tfs

    cur_outdir=outdir+"/rep"+str(run)
    os.system('mkdir -p '+ cur_outdir)
    np.savetxt(cur_outdir+'/tf_celltype_assignments.txt', z_c, fmt = '%d')
    np.savetxt(cur_outdir+'/cell_celltype_assignments.txt', c, fmt = '%d')

    # beta_P = [] #k2 by p matrix
    # M_P = [] # stores the number of regulators for each cell type
    # for i in range(k2):
    #     tf_assigned_to_celltype_i = np.where(np.array(z_c)==i)[0]
    #     Ni = np.shape(tf_assigned_to_celltype_i)[0]
    #     M_P.append(Ni)
    #     beta_P.append([(1-epsilon_P)/Ni if j in tf_assigned_to_celltype_i else epsilon_P/(p-Ni) for j in range(p)])

    # # the 2*N_P[c[i]] enables more tfs of the cell type be sampled
    # tfs_active_in_cell = np.transpose([np.random.multinomial(2*M_P[c[i]], beta_P[c[i]]) for i in range(m)])

    # see https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html
    tfs_active_in_cell_types = [np.where(z_c==i)[0] for i in range (k2)]
    rv = truncnorm(-0.8,1,loc=1,scale=1) # a_trunc = 6, b_trunc=14
    P = [[rv.rvs() if t_c_type == c_type else 0 for c_type in c] for t_c_type in z_c] 
    P = np.round(P,decimals=2)
    E = np.round(np.matmul(A,P),decimals=2)
    gene_names=np.reshape(gene,(len(gene),1))
    E_full = np.concatenate((gene_names,E),axis=1)
    np.savetxt(cur_outdir+'/P.txt', P, fmt = '%.2f', delimiter='\t')
    np.savetxt(cur_outdir+'/E.txt', E_full, fmt = '%s', delimiter='\t')

    positive_entries = E.copy()
    positive_entries[positive_entries>0] = 1

    for r in rates:
        cur_E = np.round(E * positive_entries * np.random.binomial(1, 1-r, size=(n,m)), decimals=2)
        cur_E = np.concatenate((gene_names, cur_E), axis=1)
        np.savetxt(cur_outdir+'/E_'+str(r)+'.txt', cur_E, fmt = '%s', delimiter='\t')


## for evaluating
# correlation between inferred and real TFA
# from scipy.stats import spearmanr
# testA = np.loadtxt("results/1000_n_300_p_3000_m_20_k1_10_k2_symm_dirich_5_alpha_2_pi/rep_0/A.txt")
# testP = np.loadtxt("results/1000_n_300_p_3000_m_20_k1_10_k2_symm_dirich_5_alpha_2_pi/rep_0/P.txt")
# all_corr_r2 = [np.corrcoef(np.reshape(np.concatenate((P[i,:], testP[i,:])), (2, 3000)))[0,1] for i in range(300)] # pearson
# all_corr_sp = [spearmanr(np.reshape(np.concatenate((P[i,:], testP[i,:])), (2, 3000)), axis=1)[0] for i in range(300)] # spearman
# np.where(np.sum(testP,axis=1)==0) # filter regulators. 
# np.where(np.sum(testA,axis=0)==0) # filter regulators. But from matrix A's perspective

# for jaccard
# A_reg_w_jaccard = [[np.dot(A[:,i],A[:,j])/np.sum(1*np.logical_or(A[:,i],A[:,j])) for i in range(np.shape(A)[1])] for j in range(np.shape(A)[1])]

### for testing a generic NMF
### Assuming that we use data/1000_n_300_p_3000_m_20_k1_10_k2_symm_dirich_5_alpha_2_pi/E_0.1.txt
# from sklearn.decomposition import NMF
# E_used = np.loadtxt("data/1000_n_300_p_3000_m_20_k1_10_k2_symm_dirich_5_alpha_2_pi/E_0.1.txt", usecols=range(1,m+1))
# model = NMF(n_components=10, init='nndsvd', random_state=0)
# bland_W = model.fit_transform(E_used)
# bland_H = model.components_
# bland_E = np.matmul(bland_W,bland_H)