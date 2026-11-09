# for evaluating correlation between inferred and real TFA
import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import argparse
import os

# note that test_A_Fscore_prior_sig is wrong as it samples the wrong number of edges. But we don't include them in later analysis

## Currently the function only supports nonnegative tfa
##
def compare_adj_mat_A(args):
    # testA=np.loadtxt("inferred_matrices/nmf_with_prior/rep_0/A.txt")
    if args.inferred_P_ignore_first_row_col:
        testA = (np.loadtxt(args.inferred_A, skiprows=1, dtype=str)[:, 1:]).astype(float)
    else:
        testA = np.loadtxt(args.inferred_A) # the inferred matrix, testA=np.loadtxt("inferred_matrices/nmf_with_prior/rep_0/A.txt")

    if args.transpose_inferred_A:
        testA=testA.T

    if args.inferred_P_ignore_first_row_col:
        testP = (np.loadtxt(args.inferred_P, skiprows=1, dtype=str)[:, 1:]).astype(float)
    else:
        testP = np.loadtxt(args.inferred_P) # testP = np.loadtxt("inferred_matrices/nmf_with_prior/rep_0/P.txt")

    if args.transpose_inferred_P:
        testP=testP.T

    if not args.not_scaling_inferred_P:
        # scale P rows to have unit l1 norm
        scale_testP = np.sum(np.abs(testP),axis = 1)
        testP = np.divide(testP, scale_testP[:, np.newaxis], out=np.zeros_like(testP), where=scale_testP[:, np.newaxis]!=0)
    if not args.not_scaling_inferred_A:
        testA *= scale_testP[:, np.newaxis] # will be 0 if the tf's activity is zero everywhere


    P = np.loadtxt(args.true_P)#"inputtoTFA/0.2/k20/Anc5/rep0/P.txt")


    A_lst=np.loadtxt(args.true_A, dtype=str) # ground truth but with edge-line format. "testoutput/0.2/Anc5/nw.tsv"
    prior_A_lst=np.loadtxt(args.prior_A,dtype=str) # "inputtoTFA/0.2/k20/Anc5_prior.txt"

    # "/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt"
    tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs
    gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"

    gene_dict={gene[i]:i for i in range(len(gene))}
    reg_dict={tfs[i]:i for i in range(len(tfs))}

    num_gene=len(gene)
    num_cells=P.shape[1]
    num_tfs = P.shape[0]
    A = np.zeros((num_gene,num_tfs))
    priorA = np.zeros((num_gene,num_tfs))

    for r, t in A_lst:
        A[gene_dict[t], reg_dict[r]]=1

    for r, t, _ in prior_A_lst:
        priorA[gene_dict[t], reg_dict[r]]=_

    A_nonzero = np.count_nonzero(A)
    testA_nonzero = np.count_nonzero(testA)
    top=min(args.top, testA_nonzero) 
    print("asked top "+str(args.top)+" edges. Pick "+ str(top)+" edges as the inferred network has "+str(testA_nonzero)+ " nonzero elements")

    testA_top_thr=np.partition(testA[np.nonzero(testA)], testA_nonzero-top)[testA_nonzero-top]
    testA_top = (testA>=testA_top_thr) * testA

    test_A_TP=np.count_nonzero(A*testA_top.T)

    test_A_Fscore=test_A_TP / (test_A_TP + 0.5 * (np.count_nonzero(A) + np.count_nonzero(testA_top) - 2 * test_A_TP))

    test_A_prior_TP=np.count_nonzero(priorA*testA_top.T)

    test_A_Fscore_prior=test_A_prior_TP / (test_A_prior_TP + 0.5 * (np.count_nonzero(priorA) + np.count_nonzero(testA_top) - 2 * test_A_prior_TP))

    prior_A_ground_truth_TP=np.count_nonzero(priorA*A)

    prior_A_Fscore_A =prior_A_ground_truth_TP / (prior_A_ground_truth_TP + 0.5 * (np.count_nonzero(priorA) + np.count_nonzero(A) - 2 * prior_A_ground_truth_TP))


    # permutation test
    num_it=1000
    random_A_Fscore=[]
    for _ in range(num_it):
        # A_random=np.random.binomial(1, top/num_gene/num_tfs, size=(num_gene,num_tfs))
        A_random_1d=np.zeros((A.size,1))
        A_random_1d[np.random.choice(A.size, top, replace=False)]=1
        A_random=np.reshape(A_random_1d, A.shape)
        random_A_TP=np.count_nonzero(A*A_random)
        random_A_Fscore+=[random_A_TP / (random_A_TP + 0.5 * (np.count_nonzero(A) + np.count_nonzero(A_random) - 2 * random_A_TP))]

    test_A_Fscore_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore)/len(random_A_Fscore)
    test_A_Fscore_prior_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore_prior)/len(random_A_Fscore)

    ### below for methods that estimate nonnegative tfa
    # standardize rows of testP
    # if not args.not_scaling_inferred_P:
    #     testP=testP/np.abs(testP).sum(axis=1, keepdims=True)
    # threshold by theoretical (?) baselines
    if not args.not_thresholding_inferred_P:
        testP[testP<0.5/testP.shape[1]]=0

    all_tfs_r2 = [pearsonr(P[i,:], testP[i,:]) for i in range(num_tfs)] # pearson
    all_tfs_sp = [spearmanr(np.reshape(np.concatenate((P[i,:], testP[i,:])), (2, num_cells)), axis=1) for i in range(num_tfs)] # spearman
    all_corr_r2 = [_[0] for _ in all_tfs_r2]
    all_pval_r2 = [_[1] for _ in all_tfs_r2]
    all_corr_sp = [_[0] for _ in all_tfs_sp]
    all_pval_sp = [_[1] for _ in all_tfs_sp]

    tf_prior_out_edges=np.round(np.sum(priorA, axis=0).reshape((len(tfs),1)),decimals=1)
    tf_groundtruth_out_edges=np.round(np.sum(A, axis=0).reshape((len(tfs),1)),decimals=1)
    tf_overlap_out_edges=np.round(np.sum(A*priorA, axis=0).reshape((len(tfs),1)),decimals=1)

    tf_outmat=np.concatenate((tfs.reshape((len(tfs),1)), np.array(all_corr_r2).reshape((len(tfs),1)), np.array(all_pval_r2).reshape((len(tfs),1)), np.array(all_corr_sp).reshape((len(tfs),1)), \
                    np.array(all_pval_sp).reshape((len(tfs),1)), tf_groundtruth_out_edges, tf_prior_out_edges, tf_overlap_out_edges),axis=1)
    tf_outmat=np.concatenate((np.reshape(["TF", "Pearson", "Pearson_p", "Spearman", "Spearman_p","true_out_d","prior_out_d", "overlap_out_d"],(1, tf_outmat.shape[1])), tf_outmat), axis=0)

    summary_outmat=np.array([['fscore_w_gt', 'pval_fscore_w_gt', 'fscore_w_prior', 'pval_fscore_w_prior'],[test_A_Fscore, test_A_Fscore_sig, test_A_Fscore_prior, test_A_Fscore_prior_sig]])

    np.savetxt(outdir+"/tf_specifics.txt", tf_outmat, fmt='%s',delimiter='\t')
    np.savetxt(outdir+"/stat.txt", summary_outmat, fmt='%s',delimiter='\t')

    # np.where(np.sum(testP,axis=1)==0) # filter regulators. 
    # np.where(np.sum(testA,axis=0)==0) # filter regulators. But from matrix A's perspective

#
# !!!here we assume that args.testP is tf by cells, with row names. TODO note that A not scaled by default, and P not thresholded by default
#
def compare_edge_lst_A(args): 
    
    testA = np.loadtxt(args.inferred_A, dtype=str) 
    testA = testA[testA[:, 2].astype(float).argsort()][::-1]

    ## TODO
    # if args.inferred_P_ignore_first_row_col:
    #     testP = (np.loadtxt(args.inferred_P, skiprows=1, dtype=str)[:, 1:]).astype(float)
    # else:
    #     testP = np.loadtxt(args.inferred_P) # testP = np.loadtxt("inferred_matrices/nmf_with_prior/rep_0/P.txt")

    testP = np.loadtxt(args.inferred_P, dtype=str)

    gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
    tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs

    # reorder P by tfs
    order_dict = {val: idx for idx, val in enumerate(tfs)}
    sort_indices = np.array([order_dict[val] for val in testP[:, 0]])
    testP=testP[sort_indices]
    testP=testP[:,1:].astype(float)

    # scale P rows to have unit l1 norm
    if not args.not_scaling_inferred_P:
        scale_testP = np.sum(np.abs(testP),axis = 1)
        testP /= scale_testP[:, np.newaxis]
        # testA *= scale_testP[:, np.newaxis]

    P = np.loadtxt(args.true_P)#"inputtoTFA/0.2/k20/Anc5/rep0/P.txt")

    A_lst=np.loadtxt(args.true_A, dtype=str) # ground truth but with edge-line format. "testoutput/0.2/Anc5/nw.tsv"
    prior_A_lst=np.loadtxt(args.prior_A,dtype=str) # "inputtoTFA/0.2/k20/Anc5_prior.txt"


    # A_nonzero = len(A_lst)
    testA_nonzero = len(testA)
    top=min(args.top, testA_nonzero) 
    print("asked top "+str(args.top)+" edges. Pick "+ str(top)+" edges as the inferred network has "+str(testA_nonzero)+ " nonzero elements")

    testA_top = testA[:top,:]

    A_lst_1d_view = [a+" "+b for a, b in A_lst]
    testA_top_1d_view = [a+" "+b for a, b, _ in testA_top]
    priorA_1d_view = [a+" "+b for a, b, _ in prior_A_lst]

    test_A_TP=len(np.intersect1d(A_lst_1d_view,testA_top_1d_view))

    test_A_Fscore=test_A_TP / (test_A_TP + 0.5 * (len(A_lst_1d_view) + len(testA_top) - 2 * test_A_TP))

    test_A_prior_TP=len(np.intersect1d(priorA_1d_view,testA_top_1d_view))

    test_A_Fscore_prior=test_A_prior_TP / (test_A_prior_TP + 0.5 * (len(priorA_1d_view) + len(testA_top) - 2 * test_A_prior_TP))

    # prior_A_ground_truth_TP=np.count_nonzero(priorA*A)

    # prior_A_Fscore_A =prior_A_ground_truth_TP / (prior_A_ground_truth_TP + 0.5 * (np.count_nonzero(priorA) + np.count_nonzero(A) - 2 * prior_A_ground_truth_TP))

    # permutation test
    num_it=1000
    random_A_Fscore=[]
    dummy_ground_truth=np.array(range(len(A_lst)))
    dummy_prior=np.array(range(len(prior_A_lst)))
    for _ in range(num_it):
        # A_random=np.random.binomial(1, top/num_gene/num_tfs, size=(num_gene,num_tfs))
        A_random_1d=np.random.choice(len(gene)*len(tfs), top, replace=False)
        random_A_TP=np.sum(np.isin(A_random_1d,dummy_ground_truth))
        random_A_Fscore+=[random_A_TP / (random_A_TP + 0.5 * (len(dummy_ground_truth) + top - 2 * random_A_TP))]

    test_A_Fscore_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore)/len(random_A_Fscore)
    test_A_Fscore_prior_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore_prior)/len(random_A_Fscore)

    ##  ## As tfa is allowed to take negative values, we do not do this
    # # standardize rows of testP
    # testP=testP/testP.sum(axis=1, keepdims=True)
    # # threshold by theoretical baselines
    # testP[testP<1/testP.shape[1]]=0
    num_cells=P.shape[1]
    num_tfs = P.shape[0]

    all_tfs_r2 = [pearsonr(P[i,:], testP[i,:]) for i in range(num_tfs)] # pearson
    all_tfs_sp = [spearmanr(np.reshape(np.concatenate((P[i,:], testP[i,:])), (2, num_cells)), axis=1) for i in range(num_tfs)] # spearman
    all_corr_r2 = [_[0] for _ in all_tfs_r2]
    all_pval_r2 = [_[1] for _ in all_tfs_r2]
    all_corr_sp = [_[0] for _ in all_tfs_sp]
    all_pval_sp = [_[1] for _ in all_tfs_sp]

    tf_prior_outd_dict={key:int(val) for key, val in np.vstack(np.unique(prior_A_lst[:,0], return_counts=True)).T}
    tf_prior_out_edges=np.round([tf_prior_outd_dict[r] for r in tfs],decimals=1).reshape((len(tfs),1))

    tf_gt_outd_dict={key:int(val) for key, val in np.vstack(np.unique(A_lst[:,0], return_counts=True)).T}
    tf_groundtruth_out_edges=np.round([tf_gt_outd_dict[r] for r in tfs],decimals=1).reshape((len(tfs),1))

    tf_overlap_outd_dict={key:int(val) for key, val in np.vstack(np.unique([a.split(" ", 1)[0] for a in list(set(priorA_1d_view).intersection(A_lst_1d_view))], return_counts=True)).T}
    tf_overlap_out_edges=np.round([tf_overlap_outd_dict[r] for r in tfs],decimals=1).reshape((len(tfs),1))

    tf_outmat=np.concatenate((tfs.reshape((len(tfs),1)), np.array(all_corr_r2).reshape((len(tfs),1)), np.array(all_pval_r2).reshape((len(tfs),1)), np.array(all_corr_sp).reshape((len(tfs),1)), \
                    np.array(all_pval_sp).reshape((len(tfs),1)), tf_groundtruth_out_edges, tf_prior_out_edges, tf_overlap_out_edges),axis=1)
    tf_outmat=np.concatenate((np.reshape(["TF", "Pearson", "Pearson_p", "Spearman", "Spearman_p","true_out_d","prior_out_d", "overlap_out_d"],(1, tf_outmat.shape[1])), tf_outmat), axis=0)

    summary_outmat=np.array([['fscore_w_gt', 'pval_fscore_w_gt', 'fscore_w_prior', 'pval_fscore_w_prior'],[test_A_Fscore, test_A_Fscore_sig, test_A_Fscore_prior, test_A_Fscore_prior_sig]])

    np.savetxt(outdir+"/tf_specifics.txt", tf_outmat, fmt='%s',delimiter='\t')
    np.savetxt(outdir+"/stat.txt", summary_outmat, fmt='%s',delimiter='\t')

    # np.where(np.sum(testP,axis=1)==0) # filter regulators. 
    # np.where(np.sum(testA,axis=0)==0) # filter regulators. But from matrix A's perspective




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", help="output directory. tf_specifics.txt: inferred TFA correlations to true TFA per TF and other attributes. stat.txt: other summary statistics", type=str, required=True)
    parser.add_argument("--inferred_A", help="The inferred tf by gene GRN", type=str, required=True)
    parser.add_argument("--inferred_P", help="The inferred tf by cell TFA matrix", type=str, required=True)
    parser.add_argument("--true_A", help="The true GRN edge file where each line is regulator[tab]target.", type=str, required=True)
    parser.add_argument("--prior_A", help="The prior GRN edge file used where each line is regulator[tab]target[tab]confidence.", type=str, required=True)
    parser.add_argument("--true_P", help="The true tf by cell TFA matrix", type=str, required=True)
    parser.add_argument("--tf_lst", help="a column file where each line is the name of a regulator. This specifies the order of tfs in prior network.", type=str, required=True)
    parser.add_argument("--gene", help="a column file where each line is the name of a gene. This specifies the order of genes in prior network and output expression.", type=str, required=True)
    parser.add_argument("--transpose_inferred_A", help='use this flag when the inferred GRN matrix is gene by tfs.', action='store_true')
    parser.add_argument("--transpose_inferred_P", help='use this flag when the inferred TFA matrix is cell by tfs.', action='store_true')
    parser.add_argument("--inferred_A_ignore_first_row_col", help='use this flag to ignore the first row and col of inferred A', action='store_true')
    parser.add_argument("--inferred_P_ignore_first_row_col", help='use this flag to ignore the first row and col of inferred P', action='store_true')
    parser.add_argument("--top", help="use top k edges from the inferred network", type=int, default=500)
    parser.add_argument("--A_edge_lst", help="use this flag if the inferred GRN is passed as a 2-column or 3-column (weighted) edge list for edge confidence. Default for adjacency matrix. Also, the inferred P is assumed to have row names in an unspecified order", action='store_true')
    parser.add_argument("--not_scaling_inferred_A", help='use this flag such that the inferred A is treated as it is and not depend on the inferred P. This is true for pmf-grn', action='store_true')
    parser.add_argument("--not_scaling_inferred_P", help='use this flag such that the inferred P is treated as it is and not depend on the inferred A. This is true for TIGER, which has all zero rows', action='store_true')
    parser.add_argument("--not_thresholding_inferred_P", help='use this flag such that the inferred P is not thresholded, as TFA have negative entries. This is true for regularized nca, NetREX-CF', action='store_true')

    args = parser.parse_args()
    outdir=args.outdir
    os.system("mkdir -p "+outdir)
    if args.A_edge_lst:
        compare_edge_lst_A(args)
    else:
        compare_adj_mat_A(args)


