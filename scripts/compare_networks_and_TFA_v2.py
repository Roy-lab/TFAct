# for evaluating correlation between inferred and real TFA
# !! NOTE: pval_fscore_w_prior is meant for the TFA's GRN matrix, as downstream network inference may include a larger set of tfs and genes

import numpy as np
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from scipy.stats import hypergeom
import argparse
import os

# note that test_A_Fscore_prior_sig is wrong as it samples the wrong number of edges. But we don't include them in later analysis
epsilon = 1e-14 # as TFA has scaling ambiguity, all matrices are rescaled to have unit l1 norm for each TF (depends on algorithms) and adjust the inferred GRN accordingly. This ensures that TFs with no activities will not be scaled.
## Currently the function only supports nonnegative tfa
##
def compare_adj_mat_A(args):
    print("WARNING: compare_adj_mat_A does not support evaluating both mRNA->mRNA and TFA->mRNA edge types.")
    # testA=np.loadtxt("inferred_matrices/nmf_with_prior/rep_0/A.txt")
    if args.inferred_A_ignore_first_row_col:
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
        testP = np.divide(testP, scale_testP[:, np.newaxis], out=np.zeros_like(testP), where=scale_testP[:, np.newaxis]> testP.shape[1]*epsilon)
    if not args.not_scaling_inferred_A:
        testA *= scale_testP[:, np.newaxis] # will be 0 if the tf's activity is zero everywhere

    if args.true_P:
        P = np.loadtxt(args.true_P)#"inputtoTFA/0.2/k20/Anc5/rep0/P.txt")

    A_lst=np.loadtxt(args.true_A, dtype=str) # ground truth but with edge-line format. "testoutput/0.2/Anc5/nw.tsv"
    prior_A_lst=np.loadtxt(args.prior_A,dtype=str) # "inputtoTFA/0.2/k20/Anc5_prior.txt"

    # "/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt"
    tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs
    gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
    intersected_tfs = np.intersect1d(tfs, np.unique(A_lst[:, 0]))
    intersected_genes = np.intersect1d(gene, np.unique(A_lst[:, 1]))
    intersected_loc_in_tfs = [np.array(tfs).tolist().index(t) for t in intersected_tfs]
    intersected_loc_in_genes = [np.array(gene).tolist().index(t) for t in intersected_genes]
    testA_intersected = testA[intersected_loc_in_tfs,:][:,intersected_loc_in_genes]
    int_gene_dict={intersected_genes[i]:i for i in range(len(intersected_genes))}
    int_reg_dict={intersected_tfs[i]:i for i in range(len(intersected_tfs))}

    num_cells=testP.shape[1]
    num_tfs = testP.shape[0]
    A = np.zeros((testA_intersected.shape[1], testA_intersected.shape[0]))
    priorA = np.zeros((testA_intersected.shape[1], testA_intersected.shape[0]))


    for _ in A_lst:
        r = _[0]
        t = _[1]
        if r in intersected_tfs and t in intersected_genes:
            A[int_gene_dict[t], int_reg_dict[r]]=1

    for r, t, _ in prior_A_lst:
        if r in intersected_tfs and t in intersected_genes:
            priorA[int_gene_dict[t], int_reg_dict[r]]=_

    A_nonzero = np.count_nonzero(A)
    testA_nonzero = np.count_nonzero(testA_intersected)
    top=min(args.top, testA_nonzero) 
    print("asked top "+str(args.top)+" edges. Pick "+ str(top)+" edges as the inferred network has "+str(testA_nonzero)+ " nonzero elements")

    # for evaluation on the real dataset, we need to account for negative edges, so select top edges by absolute values
    testA_top_thr=np.partition(np.abs(testA_intersected[np.nonzero(testA_intersected)]), testA_nonzero-top)[testA_nonzero-top]
    testA_top = (np.abs(testA_intersected)>=testA_top_thr) * testA_intersected

    test_A_TP=np.count_nonzero(A*testA_top.T)

    test_A_Fscore=test_A_TP / (test_A_TP + 0.5 * (np.count_nonzero(A) + np.count_nonzero(testA_top) - 2 * test_A_TP))

    test_A_prior_TP=np.count_nonzero(priorA*testA_top.T)

    test_A_Fscore_prior=test_A_prior_TP / (test_A_prior_TP + 0.5 * (np.count_nonzero(priorA) + np.count_nonzero(testA_top) - 2 * test_A_prior_TP))

    prior_A_ground_truth_TP=np.count_nonzero(priorA*A)

    prior_A_Fscore_A =prior_A_ground_truth_TP / (prior_A_ground_truth_TP + 0.5 * (np.count_nonzero(priorA) + np.count_nonzero(A) - 2 * prior_A_ground_truth_TP))


    # # permutation test
    # num_it=1000
    # random_A_Fscore=[]
    # for _ in range(num_it):
    #     # A_random=np.random.binomial(1, top/num_gene/num_tfs, size=(num_gene,num_tfs))
    #     A_random_1d=np.zeros((A.size,1))
    #     A_random_1d[np.random.choice(A.size, top, replace=False)]=1
    #     A_random=np.reshape(A_random_1d, A.shape)
    #     random_A_TP=np.count_nonzero(A*A_random)
    #     random_A_Fscore+=[random_A_TP / (random_A_TP + 0.5 * (np.count_nonzero(A) + np.count_nonzero(A_random) - 2 * random_A_TP))]

    # test_A_Fscore_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore)/len(random_A_Fscore)
    # test_A_Fscore_prior_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore_prior)/len(random_A_Fscore)
    test_A_Fscore_sig=hypergeom.sf(test_A_TP - 1, A.size, np.count_nonzero(A), top)
    test_A_Fscore_prior_sig=hypergeom.sf(test_A_prior_TP - 1, A.size, np.count_nonzero(A), top)

    ### below for methods that estimate nonnegative tfa
    # standardize rows of testP
    # if not args.not_scaling_inferred_P:
    #     testP=testP/np.abs(testP).sum(axis=1, keepdims=True)
    # threshold by theoretical (?) baselines
    if not args.not_thresholding_inferred_P:
        testP[testP<0.5/testP.shape[1]]=0

    if args.true_P:
        all_tfs_r2 = [pearsonr(P[i,:], testP[i,:]) for i in range(num_tfs)] # pearson
        all_tfs_sp = [spearmanr(np.reshape(np.concatenate((P[i,:], testP[i,:])), (2, num_cells)), axis=1) for i in range(num_tfs)] # spearman
        all_corr_r2 = [_[0] for _ in all_tfs_r2]
        all_pval_r2 = [_[1] for _ in all_tfs_r2]
        all_corr_sp = [_[0] for _ in all_tfs_sp]
        all_pval_sp = [_[1] for _ in all_tfs_sp]

        # tf_prior_out_edges=np.round(np.sum(priorA, axis=0).reshape((len(tfs),1)),decimals=1)
        # tf_groundtruth_out_edges=np.round(np.sum(A, axis=0).reshape((len(tfs),1)),decimals=1)
        # tf_overlap_out_edges=np.round(np.sum(A*priorA, axis=0).reshape((len(tfs),1)),decimals=1)

        tf_outmat=np.concatenate((tfs.reshape((len(tfs),1)), np.array(all_corr_r2).reshape((len(tfs),1)), np.array(all_pval_r2).reshape((len(tfs),1)), np.array(all_corr_sp).reshape((len(tfs),1)), \
                        np.array(all_pval_sp).reshape((len(tfs),1))),axis=1)
        tf_outmat=np.concatenate((np.reshape(["TF", "Pearson", "Pearson_p", "Spearman", "Spearman_p"],(1, tf_outmat.shape[1])), tf_outmat), axis=0)
        np.savetxt(outdir+"/tf_specifics.txt", tf_outmat, fmt='%s',delimiter='\t')

    summary_outmat=np.array([['fscore_w_gt', 'pval_fscore_w_gt', 'fscore_w_prior', 'pval_fscore_w_prior'],[test_A_Fscore, test_A_Fscore_sig, test_A_Fscore_prior, test_A_Fscore_prior_sig]])

    np.savetxt(outdir+"/stat.txt", summary_outmat, fmt='%s',delimiter='\t')

    # np.where(np.sum(testP,axis=1)==0) # filter regulators. a
    # np.where(np.sum(testA,axis=0)==0) # filter regulators. But from matrix A's perspective

def compare_edge_lst_A_edge_types(edgetype, testA, A_lst, prior_A_lst, args):
    gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
    tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs
    # A_nonzero = len(A_lst)

    intersected_tfs = np.intersect1d(tfs, np.unique(A_lst[:, 0]))
    intersected_genes = np.intersect1d(gene, np.unique(A_lst[:, 1]))
   
    A_lst = np.array([_[:2] for _ in A_lst if _[0] in intersected_tfs and _[1] in intersected_genes])

    if edgetype == 'TFA_only':
        testA = np.array([[r.removesuffix('_TFA'), t, c] for r,t,c in testA if r.removesuffix('_TFA') in intersected_tfs and r.endswith('_TFA') and t in intersected_genes and not t.endswith('_TFA')])
    elif edgetype == 'gene_only':
        testA = np.array([[r.removesuffix('_TFA'), t, c] for r,t,c in testA if r in intersected_tfs and not r.endswith('_TFA') and t in intersected_genes and not t.endswith('_TFA')])
    elif edgetype == 'max_conf':
        testA = np.array([[r.removesuffix('_TFA'), t, c] for r,t,c in testA if r.removesuffix('_TFA') in intersected_tfs and t in intersected_genes and not t.endswith('_TFA')])
        testA_tmp = dict()
        for r,s,t in testA:
            if (r,s) in testA_tmp.keys():
                if float(t) > testA_tmp[(r,s)]:
                    testA_tmp[(r,s)] = float(t)
            else:
                testA_tmp[(r,s)] = float(t)
        testA = np.array([[k[0], k[1],v] for k, v in testA_tmp.items()])
        testA = testA[testA[:, 2].astype(float).argsort()][::-1]
    else:
        raise ValueError("only support TFA_only, gene_only, and max_conf edge types.")
    
    testA_nonzero = len(testA)
    top=min(args.top, testA_nonzero) 
    print(edgetype +": asked top "+str(args.top)+" edges. Pick "+ str(top)+" for evaluation.")
    testA_top = testA[:top,:]

    A_lst_1d_view = [a+" "+b for a, b in A_lst]
    testA_top_1d_view = [a+" "+b for a, b, _ in testA_top]
    priorA_1d_view = [a+" "+b for a, b, _ in prior_A_lst]

    test_A_TP=len(np.intersect1d(A_lst_1d_view,testA_top_1d_view))

    test_A_Fscore=test_A_TP / (test_A_TP + 0.5 * (len(A_lst_1d_view) + len(testA_top) - 2 * test_A_TP))

    test_A_prior_TP=len(np.intersect1d(priorA_1d_view,testA_top_1d_view))

    test_A_Fscore_prior=test_A_prior_TP / (test_A_prior_TP + 0.5 * (len(priorA_1d_view) + len(testA_top) - 2 * test_A_prior_TP))

    test_A_Fscore_sig=hypergeom.sf(test_A_TP - 1, len(intersected_genes)*len(intersected_tfs), len(A_lst), len(testA_top))
    test_A_Fscore_prior_sig=hypergeom.sf(test_A_prior_TP - 1, len(intersected_genes)*len(intersected_tfs), len(A_lst), len(testA_top))

    summary_outmat=np.array([['fscore_w_gt', 'pval_fscore_w_gt', 'fscore_w_prior', 'pval_fscore_w_prior'],[test_A_Fscore, test_A_Fscore_sig, test_A_Fscore_prior, test_A_Fscore_prior_sig]])
    np.savetxt(outdir+"/"+ edgetype +"_" + "stat.txt", summary_outmat, fmt='%s',delimiter='\t')


#
# !!!here we assume that args.testP is tf by cells, with row names. Note that A not scaled by default, and P not thresholded by default
#
def compare_edge_lst_A(args): 
    
    testA = np.loadtxt(args.inferred_A, dtype=str) 
    testA = testA[testA[:, 2].astype(float).argsort()][::-1]
    A_lst=np.loadtxt(args.true_A, dtype=str) # ground truth but with edge-line format. "testoutput/0.2/Anc5/nw.tsv"
    prior_A_lst=np.loadtxt(args.prior_A,dtype=str) # "inputtoTFA/0.2/k20/Anc5_prior.txt"

    ## TODO
    # if args.inferred_P_ignore_first_row_col:
    #     testP = (np.loadtxt(args.inferred_P, skiprows=1, dtype=str)[:, 1:]).astype(float)
    # else:
    #     testP = np.loadtxt(args.inferred_P) # testP = np.loadtxt("inferred_matrices/nmf_with_prior/rep_0/P.txt")

    if args.separate_eval_by_edge_types: # have nothing to do with P
        compare_edge_lst_A_edge_types("gene_only", testA, A_lst, prior_A_lst, args)
        compare_edge_lst_A_edge_types("TFA_only", testA, A_lst, prior_A_lst, args)
        compare_edge_lst_A_edge_types("max_conf", testA, A_lst, prior_A_lst, args)
    
    else:
        gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
        tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs

        if args.inferred_P is not None:
            testP = np.loadtxt(args.inferred_P, dtype=str)
            # reorder P by tfs
            order_dict = {val: idx for idx, val in enumerate(tfs)}
            # activity of not active tf are set to 0
            non_active_TF = np.concatenate((np.setdiff1d(tfs,testP[:,0]).reshape((-1,1)), np.array([['0' for i in range(testP.shape[1]-1)] for j in range(len(np.setdiff1d(tfs,testP[:,0])))])), axis=1)
            testP = np.concatenate((testP, non_active_TF), axis=0)

            sort_indices = np.array([order_dict[val] for val in testP[:, 0]])
            testP=testP[sort_indices]
            testP=testP[:,1:].astype(float)

            # scale P rows to have unit l1 norm
            if not args.not_scaling_inferred_P:
                scale_testP = np.sum(np.abs(testP),axis = 1)
                testP = np.divide(testP, scale_testP[:, np.newaxis], out=np.zeros_like(testP), where=scale_testP[:, np.newaxis]> testP.shape[1]*epsilon)
                # testA *= scale_testP[:, np.newaxis]

        if args.true_P:
            P = np.loadtxt(args.true_P)#"inputtoTFA/0.2/k20/Anc5/rep0/P.txt")

        # A_nonzero = len(A_lst)

        intersected_tfs = np.intersect1d(tfs, np.unique(A_lst[:, 0]))
        intersected_genes = np.intersect1d(gene, np.unique(A_lst[:, 1]))
        # prior_intersected_tfs = np.intersect1d(tfs, np.unique(prior_A_lst[:, 0]))
        # prior_intersected_genes = np.intersect1d(gene, np.unique(prior_A_lst[:, 1]))

        A_lst = np.array([_[:2] for _ in A_lst if _[0] in intersected_tfs and _[1] in intersected_genes])

        # testA_top_for_prior = testA[:top,:]
        testA = np.array([[r, t, c] for r,t,c in testA if r in intersected_tfs and t in intersected_genes])

        testA_nonzero = len(testA)
        top=min(args.top, testA_nonzero) 
        print("asked top "+str(args.top)+" edges. Pick "+ str(top)+" for evaluation.")
        testA_top = testA[:top,:]

        A_lst_1d_view = [a+" "+b for a, b in A_lst]
        testA_top_1d_view = [a+" "+b for a, b, _ in testA_top]
        # testA_top_1d_view_for_prior = [a+" "+b for a, b, _ in testA_top_for_prior] # not used for this. We focus on priors recovered in terms of the gs
        priorA_1d_view = [a+" "+b for a, b, _ in prior_A_lst]

        test_A_TP=len(np.intersect1d(A_lst_1d_view,testA_top_1d_view))

        test_A_Fscore=test_A_TP / (test_A_TP + 0.5 * (len(A_lst_1d_view) + len(testA_top) - 2 * test_A_TP))

        test_A_prior_TP=len(np.intersect1d(priorA_1d_view,testA_top_1d_view))

        test_A_Fscore_prior=test_A_prior_TP / (test_A_prior_TP + 0.5 * (len(priorA_1d_view) + len(testA_top) - 2 * test_A_prior_TP))

        # prior_A_ground_truth_TP=np.count_nonzero(priorA*A)

        # prior_A_Fscore_A =prior_A_ground_truth_TP / (prior_A_ground_truth_TP + 0.5 * (np.count_nonzero(priorA) + np.count_nonzero(A) - 2 * prior_A_ground_truth_TP))

        # # permutation test
        # num_it=1000
        # random_A_Fscore=[]
        # dummy_ground_truth=np.array(range(len(A_lst)))
        # top = len(testA_top)
        # #dummy_prior=np.array(range(len(prior_A_lst)))
        # for _ in range(num_it):
        #     # A_random=np.random.binomial(1, top/num_gene/num_tfs, size=(num_gene,num_tfs))
        #     A_random_1d=np.random.choice(len(intersected_genes)*len(intersected_tfs), top, replace=False)
        #     random_A_TP=np.sum(np.isin(A_random_1d,dummy_ground_truth))
        #     random_A_Fscore+=[random_A_TP / (random_A_TP + 0.5 * (len(dummy_ground_truth) + top - 2 * random_A_TP))]
        # test_A_Fscore_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore)/len(random_A_Fscore)
        # test_A_Fscore_prior_sig=sum(np.array(random_A_Fscore) >= test_A_Fscore_prior)/len(random_A_Fscore)
        test_A_Fscore_sig=hypergeom.sf(test_A_TP - 1, len(intersected_genes)*len(intersected_tfs), len(A_lst), len(testA_top))
        test_A_Fscore_prior_sig=hypergeom.sf(test_A_prior_TP - 1, len(intersected_genes)*len(intersected_tfs), len(prior_A_lst), len(testA_top))

        ##  ## As tfa is allowed to take negative values, we do not do this
        # # standardize rows of testP
        # testP=testP/testP.sum(axis=1, keepdims=True)
        # # threshold by theoretical baselines
        # testP[testP<1/testP.shape[1]]=0
        if args.true_P:
            num_cells=P.shape[1]
            num_tfs = P.shape[0]

            all_tfs_r2 = [pearsonr(P[i,:], testP[i,:]) for i in range(num_tfs)] # pearson
            all_tfs_sp = [spearmanr(np.reshape(np.concatenate((P[i,:], testP[i,:])), (2, num_cells)), axis=1) for i in range(num_tfs)] # spearman
            all_corr_r2 = [_[0] for _ in all_tfs_r2]
            all_pval_r2 = [_[1] for _ in all_tfs_r2]
            all_corr_sp = [_[0] for _ in all_tfs_sp]
            all_pval_sp = [_[1] for _ in all_tfs_sp]
            # tf_prior_outd_dict={key:int(val) for key, val in np.vstack(np.unique(prior_A_lst[:,0], return_counts=True)).T}
            # tf_prior_out_edges=np.round([tf_prior_outd_dict[r] for r in tfs],decimals=1).reshape((len(tfs),1))

            # tf_gt_outd_dict={key:int(val) for key, val in np.vstack(np.unique(A_lst[:,0], return_counts=True)).T}
            # tf_groundtruth_out_edges=np.round([tf_gt_outd_dict[r] for r in tfs],decimals=1).reshape((len(tfs),1))

            # tf_overlap_outd_dict={key:int(val) for key, val in np.vstack(np.unique([a.split(" ", 1)[0] for a in list(set(priorA_1d_view).intersection(A_lst_1d_view))], return_counts=True)).T}
            # tf_overlap_out_edges=np.round([tf_overlap_outd_dict[r] for r in tfs],decimals=1).reshape((len(tfs),1))

            tf_outmat=np.concatenate((tfs.reshape((len(tfs),1)), np.array(all_corr_r2).reshape((len(tfs),1)), np.array(all_pval_r2).reshape((len(tfs),1)), np.array(all_corr_sp).reshape((len(tfs),1)), \
                            np.array(all_pval_sp).reshape((len(tfs),1))),axis=1)
            tf_outmat=np.concatenate((np.reshape(["TF", "Pearson", "Pearson_p", "Spearman", "Spearman_p"],(1, tf_outmat.shape[1])), tf_outmat), axis=0)
            np.savetxt(outdir+"/tf_specifics.txt", tf_outmat, fmt='%s',delimiter='\t')

        summary_outmat=np.array([['fscore_w_gt', 'pval_fscore_w_gt', 'fscore_w_prior', 'pval_fscore_w_prior'],[test_A_Fscore, test_A_Fscore_sig, test_A_Fscore_prior, test_A_Fscore_prior_sig]])

        
        np.savetxt(outdir+"/stat.txt", summary_outmat, fmt='%s',delimiter='\t')

    # np.where(np.sum(testP,axis=1)==0) # filter regulators. 
    # np.where(np.sum(testA,axis=0)==0) # filter regulators. But from matrix A's perspective




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", help="output directory. tf_specifics.txt: inferred TFA correlations to true TFA per TF and other attributes. stat.txt: other summary statistics", type=str, required=True)
    parser.add_argument("--inferred_A", help="The inferred tf by gene GRN", type=str, required=True)
    parser.add_argument("--inferred_P", help="The inferred tf by cell TFA matrix", type=str)
    parser.add_argument("--true_A", help="The true GRN edge file where each line is regulator[tab]target.", type=str, required=True)
    parser.add_argument("--prior_A", help="The prior GRN edge file used where each line is regulator[tab]target[tab]confidence.", type=str, required=True)
    parser.add_argument("--true_P", help="The true tf by cell TFA matrix", type=str)
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

    parser.add_argument("--separate_eval_by_edge_types", help="use this flag if the inferred GRN has TFA->mRNA edges. We will evaluate: 1. TFA->mRNA only, 2. mRNA->mRNA, and 3. max confidence of 1 and 2. Other edge types (e.g. mRNA->TFA) will be discarded.", action='store_true')

    args = parser.parse_args()
    outdir=args.outdir
    os.system("mkdir -p "+outdir)
    if args.A_edge_lst:
        compare_edge_lst_A(args)
    else:
        compare_adj_mat_A(args)


