import numpy as np
import pandas as pd
import argparse

def reformat_tfa(args):
    exp = pd.read_csv(args.exp, sep='\t', index_col=0)
    tfa = pd.read_csv(args.tfa, sep=args.sep, header=None)
    tf_lst = pd.read_csv(args.tf_lst, header=None)

    tfa.index = tf_lst.iloc[:,0]
    tfa.columns = [ "C{}".format(i) for i in range(tfa.shape[1])]

    tfa_filtered = tfa.loc[:, exp.columns]
    tfa_filtered.T.to_csv(args.out_tfa, sep='\t', index=True, header=True, compression="gzip", index_label="")

def reformat_tfa_first_col_is_TF(args):
    exp = pd.read_csv(args.exp, sep='\t', index_col=0)
    tfa = pd.read_csv(args.tfa, sep=args.sep, header=None)

    tfa.index = tfa.iloc[:,0]
    tfa = tfa.iloc[:,1:]
    tfa.columns = [ "C{}".format(i) for i in range(tfa.shape[1])]

    tfa_filtered = tfa.loc[:, exp.columns]
    tfa_filtered.T.to_csv(args.out_tfa, sep='\t', index=True, header=True, compression="gzip", index_label="")

def reformat_tfa_first_col_is_TF_fill_missing_TFA_w_exp(args):
    exp = pd.read_csv(args.exp, sep='\t', index_col=0)
    tfa = pd.read_csv(args.tfa, sep=args.sep, header=None)
    tf_lst = pd.read_csv(args.tf_lst, header=None)

    tfa.index = tfa.iloc[:,0]
    tfa = tfa.iloc[:,1:]
    tfa.columns = [ "C{}".format(i) for i in range(tfa.shape[1])]

    tf_wo_tfa = np.intersect1d(np.setdiff1d(tf_lst.iloc[:,0], tfa.index),exp.index)

    tfa_filtered = tfa.loc[:, exp.columns]
    tf_wo_tfa_exp = exp.loc[tf_wo_tfa,:]
    tfa_out = pd.concat([tfa_filtered, tf_wo_tfa_exp]).sort_index()
    tfa_out.T.to_csv(args.out_tfa, sep='\t', index=True, header=True, compression="gzip", index_label="")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp", help="a tsv expression matrix. Headers correspond to cell names", type=str, required=True)
    parser.add_argument("--tf_lst", help="a column file where each line is the name of a regulator. If not supplied, the first col of exp is TFs.", type=str)
    parser.add_argument("--algo_name", help="Name of the tfa inference algorithm.", type=str)
    parser.add_argument("--tfa", help="tfa file, without a header and assumes different formats for different algos. The missing header is assumed to be C0\\tC1...", type=str, required=True)
    parser.add_argument("--out_tfa", help="output file, sample by tsv", type=str, required=True)
    parser.add_argument("--sep", help="tfa file's separator. Default ' '. Pass $'\\t' for tsv", type=str, default=' ')
    parser.add_argument("--replace_missing_TFA_w_exp", help="add this flag if you want to treat missing TFA as expression", action='store_true')
    args = parser.parse_args()
    if args.tf_lst is not None:
        if args.replace_missing_TFA_w_exp:
            reformat_tfa_first_col_is_TF_fill_missing_TFA_w_exp(args)
        else:
            reformat_tfa(args)
    else:
        reformat_tfa_first_col_is_TF(args)
