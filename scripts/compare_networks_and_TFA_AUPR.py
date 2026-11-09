# for evaluating correlation between inferred and real TFA
import numpy as np
# from scipy.stats import spearmanr
# from scipy.stats import pearsonr
from sklearn.metrics import precision_recall_curve, auc
import argparse
import os

# note that test_A_Fscore_prior_sig is wrong as it samples the wrong number of edges. But we don't include them in later analysis
epsilon = 1e-14 # as TFA has scaling ambiguity, all matrices are rescaled to have unit l1 norm for each TF (depends on algorithms) and adjust the inferred GRN accordingly. This ensures that TFs with no activities will not be scaled.
## Currently the function only supports nonnegative tfa
##
def compare_adj_mat_A(args):
    print("WARNING: compare_adj_mat_A does not support evaluating only mRNA->mRNA or only TFA->mRNA edge types.")
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


    A_lst = np.loadtxt(args.true_A, dtype=str, usecols=(0, 1)) # ground truth but with edge-line format. "testoutput/0.2/Anc5/nw.tsv"

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
    int_gene_id_dict={i:intersected_genes[i] for i in range(len(intersected_genes))}
    int_reg_id_dict={i:intersected_tfs[i] for i in range(len(intersected_tfs))}

    A = np.zeros((len(intersected_tfs), len(intersected_genes)))

    for _ in A_lst:
        r = _[0]
        t = _[1]
        if r in intersected_tfs and t in intersected_genes:
            A[int_reg_dict[r], int_gene_dict[t]]=1

    # regulator-specific AUPR
    gt_reg_cnt = np.array(np.unique(A_lst[:,0], return_counts=True)).T
    gt_reg_cnt = gt_reg_cnt[np.where(np.array(gt_reg_cnt[:,1], dtype=int) > args.tf_specific_AUPR_filter_ntarget_thr)]
    gt_reg_cnt = gt_reg_cnt[np.where(np.isin(gt_reg_cnt[:,0],intersected_tfs))]
    gt_reg_filtered = gt_reg_cnt[:,0]
    reg_specific_aupr = []

    for r in gt_reg_filtered:
        precision, recall, _ = precision_recall_curve(A[int_reg_dict[r], :], np.abs(testA_intersected[int_reg_dict[r], :]))
        reg_specific_aupr.append(auc(recall, precision))
    
    with open(outdir + "/AUPR_TF.txt", 'w') as f:
        for i in range(len(reg_specific_aupr)):
            f.write(gt_reg_filtered[i] + '\t' + str(reg_specific_aupr[i])+ '\n')

    # overall aupr
    A = A.reshape(A.shape[0] * A.shape[1])
    inferred_A = testA_intersected.reshape(testA_intersected.shape[0] * testA_intersected.shape[1])

    precision, recall, _ = precision_recall_curve(A, np.abs(inferred_A))
    aupr = auc(recall, precision)

    with open(outdir + "/AUPR.txt", 'w') as f:
        f.write(str(aupr))
    
    if not args.not_outputing_top_edges:
        testA_nonzero = np.count_nonzero(testA_intersected)
        top=min(args.top, testA_nonzero) 
        print("asked top "+str(args.top)+" edges. Pick "+ str(top)+" edges as the filtered inferred network has "+str(testA_nonzero)+ " nonzero elements")

        # for evaluation on the real dataset, we need to account for negative edges, so select top edges by absolute values
        testA_top_thr=np.partition(np.abs(testA_intersected[np.nonzero(testA_intersected)]), testA_nonzero-top)[testA_nonzero-top]
        testA_top = (np.abs(testA_intersected)>=testA_top_thr) * testA_intersected

        out_A = outdir+"/A_inferred_top_" + str(args.top) + ".txt"
        with open(out_A, 'w') as f:
            for r,t in zip(np.nonzero(testA_top)[0], np.nonzero(testA_top)[1]):
                f.write(int_reg_id_dict[r] +"\t"+ int_gene_id_dict[t]+"\n")
        
        with open(outdir + "/Enr_filelist_top_"+ str(args.top) + ".txt", 'w') as f:
            f.write("gt\t"+args.true_A+"\t1\n")
            f.write("inferred\t"+out_A+"\t2")


#
# !!!here we assume that args.testP is tf by cells, with row names. Note that A not scaled by default, and P not thresholded by default
#
def compare_edge_lst_A(args): 
    
    testA = np.loadtxt(args.inferred_A, dtype=str) 
    testA = testA[testA[:, 2].astype(float).argsort()][::-1]

    testP = np.loadtxt(args.inferred_P, dtype=str)

    gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
    tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs

    # reorder P by tfs
    order_dict = {val: idx for idx, val in enumerate(tfs)}

    # # activity of not active tf are set to 0
    # non_active_TF = np.concatenate((np.setdiff1d(tfs,testP[:,0]).reshape((-1,1)), np.array([['0' for i in range(testP.shape[1]-1)] for j in range(len(np.setdiff1d(tfs,testP[:,0])))])), axis=1)
    # testP = np.concatenate((testP, non_active_TF), axis=0)
    missing_tfs = np.setdiff1d(tfs, testP[:, 0])

    non_active_TF = np.column_stack((missing_tfs, np.full((len(missing_tfs), testP.shape[1] - 1), '0')))

    testP = np.concatenate((testP, non_active_TF), axis=0)

    sort_indices = np.array([order_dict[val] for val in testP[:, 0]])
    testP=testP[sort_indices]
    testP=testP[:,1:].astype(float)

    # scale P rows to have unit l1 norm
    if not args.not_scaling_inferred_P:
        scale_testP = np.sum(np.abs(testP),axis = 1)
        testP = np.divide(testP, scale_testP[:, np.newaxis], out=np.zeros_like(testP), where=scale_testP[:, np.newaxis]> testP.shape[1]*epsilon)

    A_lst=np.loadtxt(args.true_A, dtype=str, usecols=(0, 1)) # ground truth but with edge-line format. "testoutput/0.2/Anc5/nw.tsv"    

    # A_nonzero = len(A_lst)
    
    intersected_tfs = np.intersect1d(tfs, np.unique(A_lst[:, 0]))
    intersected_genes = np.intersect1d(gene, np.unique(A_lst[:, 1]))

    testA = np.array([[r, t, c] for r,t,c in testA if r in intersected_tfs and t in intersected_genes])

    A_lst = np.array([_[:2] for _ in A_lst if _[0] in intersected_tfs and _[1] in intersected_genes])

    # AUPR, in terms of edges formed by ground truth's tfs (intersected with tfs provided as input to algos) and genes (intersected with genes provided as input to algos)
    tf_dict = {intersected_tfs[i]:i for i in range(len(intersected_tfs))}
    gene_dict = {intersected_genes[i]:i for i in range(len(intersected_genes))}

    A = np.zeros((len(intersected_tfs), len(intersected_genes)))
    inferred_A = np.zeros((len(intersected_tfs), len(intersected_genes)))
    for _ in A_lst:
        r = _[0]
        t = _[1]
        A[tf_dict[r], gene_dict[t]]=1

    for r, t, c in testA:
        inferred_A[tf_dict[r], gene_dict[t]]=c
    
    # regulator-specific AUPR
    gt_reg_cnt = np.array(np.unique(A_lst[:,0], return_counts=True)).T
    gt_reg_cnt = gt_reg_cnt[np.where(np.array(gt_reg_cnt[:,1], dtype=int) > args.tf_specific_AUPR_filter_ntarget_thr)]
    gt_reg_cnt = gt_reg_cnt[np.where(np.isin(gt_reg_cnt[:,0],intersected_tfs))]
    gt_reg_filtered = gt_reg_cnt[:,0]
    reg_specific_aupr = []

    for r in gt_reg_filtered:
        precision, recall, _ = precision_recall_curve(A[tf_dict[r], :], np.abs(inferred_A[tf_dict[r], :]))
        reg_specific_aupr.append(auc(recall, precision))
    
    with open(outdir + "/AUPR_TF.txt", 'w') as f:
        for i in range(len(reg_specific_aupr)):
            f.write(gt_reg_filtered[i] + '\t' + str(reg_specific_aupr[i])+'\n')


    # all AUPR
    A = A.reshape(A.shape[0] * A.shape[1])
    inferred_A = inferred_A.reshape(inferred_A.shape[0] * inferred_A.shape[1])

    precision, recall, _ = precision_recall_curve(A, np.abs(inferred_A))
    aupr = auc(recall, precision)

    with open(outdir + "/AUPR.txt", 'w') as f:
        f.write(str(aupr))

    if not args.not_outputing_top_edges:
        testA_nonzero = len(testA)
        top=min(args.top, testA_nonzero) 
        print("asked top "+str(args.top)+" edges. Pick "+ str(top)+" edges as the filtered inferred network has "+str(testA_nonzero)+ " nonzero elements")

        testA_top = testA[:top,:]
        out_A = outdir+"/A_inferred_top_" + str(args.top) + ".txt"
        np.savetxt(out_A, testA_top, fmt='%s',delimiter='\t')
        with open(outdir + "/Enr_filelist_top_"+ str(args.top) + ".txt", 'w') as f:
            f.write("gt\t"+args.true_A+"\t1\n")
            f.write("inferred\t"+out_A+"\t2")


def compare_edge_lst_not_TFA_A_edge_types(edgetype, testA, A_lst, args):
    gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
    tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs

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
    
    # AUPR, in terms of edges formed by ground truth's tfs (intersected with tfs provided as input to algos) and genes (intersected with genes provided as input to algos)
    tf_dict = {intersected_tfs[i]:i for i in range(len(intersected_tfs))}
    gene_dict = {intersected_genes[i]:i for i in range(len(intersected_genes))}

    A = np.zeros((len(intersected_tfs), len(intersected_genes)))
    inferred_A = np.zeros((len(intersected_tfs), len(intersected_genes)))
    for _ in A_lst:
        r = _[0]
        t = _[1]
        A[tf_dict[r], gene_dict[t]]=1

    for r, t, c in testA:
        inferred_A[tf_dict[r], gene_dict[t]]=float(c)
    
    # regulator-specific AUPR: include regulators that have a sufficient number of targets only
    gt_reg_cnt = np.array(np.unique(A_lst[:,0], return_counts=True)).T
    gt_reg_cnt = gt_reg_cnt[np.where(np.array(gt_reg_cnt[:,1], dtype=int) > args.tf_specific_AUPR_filter_ntarget_thr)]
    gt_reg_cnt = gt_reg_cnt[np.where(np.isin(gt_reg_cnt[:,0],intersected_tfs))]
    gt_reg_filtered = gt_reg_cnt[:,0]
    reg_specific_aupr = []

    for r in gt_reg_filtered:
        precision, recall, _ = precision_recall_curve(A[tf_dict[r], :], np.abs(inferred_A[tf_dict[r], :]))
        reg_specific_aupr.append(auc(recall, precision))
    
    with open(outdir + "/"+ edgetype +"_" +"AUPR_TF.txt", 'w') as f:
        for i in range(len(reg_specific_aupr)):
            f.write(gt_reg_filtered[i] + '\t' + str(reg_specific_aupr[i])+'\n')


    # all AUPR
    A = A.reshape(A.shape[0] * A.shape[1])
    inferred_A = inferred_A.reshape(inferred_A.shape[0] * inferred_A.shape[1])

    precision, recall, _ = precision_recall_curve(A, np.abs(inferred_A))
    aupr = auc(recall, precision)

    with open(outdir + "/"+ edgetype +"_" +"AUPR.txt", 'w') as f:
        f.write(str(aupr))

    
    if not args.not_outputing_top_edges:
        testA_nonzero = len(testA)
        top=min(args.top, testA_nonzero) 
        print(edgetype +": asked top "+str(args.top)+" edges. Pick "+ str(top)+" edges as the filtered inferred network has "+str(testA_nonzero)+ " nonzero elements")

        testA_top = testA[:top,:]
        out_A = outdir+ "/"+ edgetype +"_" +"A_inferred_top_" + str(args.top) + ".txt"
        np.savetxt(out_A, testA_top, fmt='%s',delimiter='\t')
        with open(outdir + "/"+ edgetype +"_" + "Enr_filelist_top_"+ str(args.top) + ".txt", 'w') as f:
            f.write("gt\t"+args.true_A+"\t1\n")
            f.write("inferred\t"+out_A+"\t2")
    

def compare_edge_lst_A_not_TFA(args): 
    testA = np.loadtxt(args.inferred_A, dtype=str) 
    testA = testA[testA[:, 2].astype(float).argsort()][::-1]
    A_lst=np.loadtxt(args.true_A, dtype=str, usecols=(0, 1)) # ground truth but with edge-line format. "testoutput/0.2/Anc5/nw.tsv"

    if args.separate_eval_by_edge_types:
        compare_edge_lst_not_TFA_A_edge_types("gene_only", testA, A_lst, args)
        compare_edge_lst_not_TFA_A_edge_types("TFA_only", testA, A_lst, args)
        compare_edge_lst_not_TFA_A_edge_types("max_conf", testA, A_lst, args)
    else:
        gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
        tfs=np.loadtxt(args.tf_lst, dtype=str) # which requires to know the order of genes and tfs

        intersected_tfs = np.intersect1d(tfs, np.unique(A_lst[:, 0]))
        intersected_genes = np.intersect1d(gene, np.unique(A_lst[:, 1]))

        testA = np.array([[r, t, c] for r,t,c in testA if r in intersected_tfs and t in intersected_genes])

        A_lst = np.array([_[:2] for _ in A_lst if _[0] in intersected_tfs and _[1] in intersected_genes])

        # AUPR, in terms of edges formed by ground truth's tfs (intersected with tfs provided as input to algos) and genes (intersected with genes provided as input to algos)
        tf_dict = {intersected_tfs[i]:i for i in range(len(intersected_tfs))}
        gene_dict = {intersected_genes[i]:i for i in range(len(intersected_genes))}

        A = np.zeros((len(intersected_tfs), len(intersected_genes)))
        inferred_A = np.zeros((len(intersected_tfs), len(intersected_genes)))
        for _ in A_lst:
            r = _[0]
            t = _[1]
            A[tf_dict[r], gene_dict[t]]=1

        for r, t, c in testA:
            inferred_A[tf_dict[r], gene_dict[t]]=float(c)
        
        # regulator-specific AUPR: include regulators that have a sufficient number of targets only
        gt_reg_cnt = np.array(np.unique(A_lst[:,0], return_counts=True)).T
        gt_reg_cnt = gt_reg_cnt[np.where(np.array(gt_reg_cnt[:,1], dtype=int) > args.tf_specific_AUPR_filter_ntarget_thr)]
        gt_reg_cnt = gt_reg_cnt[np.where(np.isin(gt_reg_cnt[:,0],intersected_tfs))]
        gt_reg_filtered = gt_reg_cnt[:,0]
        reg_specific_aupr = []

        for r in gt_reg_filtered:
            precision, recall, _ = precision_recall_curve(A[tf_dict[r], :], np.abs(inferred_A[tf_dict[r], :]))
            reg_specific_aupr.append(auc(recall, precision))
        
        with open(outdir + "/AUPR_TF.txt", 'w') as f:
            for i in range(len(reg_specific_aupr)):
                f.write(gt_reg_filtered[i] + '\t' + str(reg_specific_aupr[i])+'\n')


        # all AUPR
        A = A.reshape(A.shape[0] * A.shape[1])
        inferred_A = inferred_A.reshape(inferred_A.shape[0] * inferred_A.shape[1])

        precision, recall, _ = precision_recall_curve(A, np.abs(inferred_A))
        aupr = auc(recall, precision)

        with open(outdir + "/AUPR.txt", 'w') as f:
            f.write(str(aupr))

        if not args.not_outputing_top_edges:
            testA_nonzero = len(testA)
            top=min(args.top, testA_nonzero) 
            print("asked top "+str(args.top)+" edges. Pick "+ str(top)+" edges as the filtered inferred network has "+str(testA_nonzero)+ " nonzero elements")

            testA_top = testA[:top,:]
            out_A = outdir+"/A_inferred_top_" + str(args.top) + ".txt"
            np.savetxt(out_A, testA_top, fmt='%s',delimiter='\t')
            with open(outdir + "/Enr_filelist_top_"+ str(args.top) + ".txt", 'w') as f:
                f.write("gt\t"+args.true_A+"\t1\n")
                f.write("inferred\t"+out_A+"\t2")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", help="output directory. AUPR.txt: a single value. A_inferred_top_$\{top\}.txt: 3 column tab delim (reg, target, score). Enr_filelist_top_$\{top\}.txt: the 2 row, 3 column (network_file, path, id)", type=str, required=True)
    parser.add_argument("--inferred_A", help="The inferred tf by gene GRN", type=str, required=True)
    parser.add_argument("--inferred_P", help="The inferred tf by cell TFA matrix", type=str)
    parser.add_argument("--true_A", help="The true GRN edge file where each line is regulator[tab]target.", type=str, required=True)
    parser.add_argument("--tf_lst", help="a column file where each line is the name of a regulator. This specifies the order of tfs in prior network.", type=str, required=True)
    parser.add_argument("--gene", help="a column file where each line is the name of a gene. This specifies the order of genes in prior network and output expression.", type=str, required=True)
    parser.add_argument("--transpose_inferred_A", help='use this flag when the inferred GRN matrix is gene by tfs.', action='store_true')
    parser.add_argument("--transpose_inferred_P", help='use this flag when the inferred TFA matrix is cell by tfs.', action='store_true')
    parser.add_argument("--inferred_A_ignore_first_row_col", help='use this flag to ignore the first row and col of inferred A', action='store_true')
    parser.add_argument("--inferred_P_ignore_first_row_col", help='use this flag to ignore the first row and col of inferred P', action='store_true')
    parser.add_argument("--top", help="use top k edges from the inferred network for outputing A_inferred_intersected.txt", type=int, default=500)
    parser.add_argument("--not_outputing_top_edges", help='use this flag to not output the top edges for predictable TFs.', action='store_true')
    parser.add_argument("--A_edge_lst", help="use this flag if the inferred GRN is passed as a 2-column or 3-column (weighted) edge list for edge confidence. Default for adjacency matrix. Also, the inferred P is assumed to have row names in an unspecified order", action='store_true')
    parser.add_argument("--not_scaling_inferred_A", help='use this flag such that the inferred A is treated as it is and not depend on the inferred P. This is true for pmf-grn', action='store_true')
    parser.add_argument("--not_scaling_inferred_P", help='use this flag such that the inferred P is treated as it is and not depend on the inferred A. This is true for TIGER, which has all zero rows', action='store_true')
    parser.add_argument("--not_thresholding_inferred_P", help='[Not implemented as P is not evaluated] use this flag such that the inferred P is not thresholded, as TFA have negative entries. This is true for regularized nca, NetREX-CF', action='store_true')
    parser.add_argument("--not_outputing_AUPR", help='use this flag such AUPR is not output. This is used for dealing snakemake issue', action='store_true')
    parser.add_argument("--tf_specific_AUPR_filter_ntarget_thr", help='For computing regulator-specific AUPR, remove TFs that have less than thr number of targets', type=int, default=10)

    parser.add_argument("--separate_eval_by_edge_types", help="use this flag if the inferred GRN has TFA->mRNA edges. We will evaluate: 1. TFA->mRNA only, 2. mRNA->mRNA, and 3. max confidence of 1 and 2. Other edge types (e.g. mRNA->TFA) will be discarded.", action='store_true')

    args = parser.parse_args()
    outdir=args.outdir
    os.system("mkdir -p "+outdir)
    if args.A_edge_lst:
        if args.inferred_P is not None:
            compare_edge_lst_A(args)
        else:
            compare_edge_lst_A_not_TFA(args)
    else:
        compare_adj_mat_A(args)


