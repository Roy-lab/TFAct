import argparse
import numpy as np
import pandas as pd

def convert_prior(args): 
    # tf_lst = pd.read_csv("/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/cluster1_allregulators.txt", header=None)
    # exp = pd.read_csv("/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/subsamples/cluster1/dataset.0.txt", index_col=0, sep='\t').T
    # # prior_A_lst = pd.read_csv("/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/cluster1_network.txt", sep='\t', header=None) 
    # tfa = pd.read_csv("/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/subsamples/cluster1/TFA/ibpg_2fct/0.tsv.gz", sep='\t', index_col=0)

    exp = pd.read_csv(args.exp, index_col=0, sep='\t').T
    # prior_A_lst = pd.read_csv(args.prior_A, sep='\t', header=None) 
    if args.expression_only:
        exp.to_csv(args.out, sep='\t', index=False, header=True)
        return
    tf_lst = pd.read_csv(args.tf_lst, header=None)
    
    if args.tfa_transposed_no_header:
        tfa = pd.read_csv(args.tfa, sep='\t', index_col=0, header = None).T
        tfa.index = ["C"+str(_) for _ in range(tfa.shape[0])]
    else:
        tfa = pd.read_csv(args.tfa, sep='\t', index_col=0)
    tfa = tfa.loc[exp.index, np.intersect1d(tf_lst[0], tfa.columns)]
    if args.replace:
        exp.update(tfa)
    else:
        tfa.columns = [_ + "_TFA" for _ in tfa.columns]
        exp = pd.concat([exp, tfa], axis=1)
    
    exp.to_csv(args.out, sep='\t', index=False, header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp", help="a tsv expression matrix (gene by sample) with row names and col names.", type=str, required=True)
    parser.add_argument("--tf_lst", help="a column file where each line is the name of a regulator.", type=str)
    # parser.add_argument("--prior_A", help="The prior GRN edge file used where each line is regulator[tab]target[tab]confidence.", type=str, required=True)
    parser.add_argument("--out", help="Path of an output tsv expression matrix (sample by gene) without row names as required by merlin.", type=str, required=True)
    parser.add_argument("--replace", help="Add this flag if want to replace the expression matrix's TF expression with TFA (only TFs in both tf_lst and tfa). Otherwise concat these TFA features.", action="store_true")
    parser.add_argument("--tfa", help="sample by tfa tsv file path", type=str)
    parser.add_argument("--tfa_transposed_no_header", help="Add this flag if processing a gene by sample TFA without a header", action="store_true")
    parser.add_argument("--expression_only", help="just reformat the expression to merlin input.", action="store_true")
    
    args = parser.parse_args()
    convert_prior(args)