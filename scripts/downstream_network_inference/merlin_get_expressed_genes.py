import argparse
import numpy as np
import pandas as pd

def get_expressed_genes(args): 
    exp = pd.read_csv(args.exp, index_col=0, sep='\t', header=None)
    idx = exp.index[(exp > 0).sum(axis=1) >= args.thr]
    with open(args.out, 'w') as f:
        for i in idx:
            _ = f.write(f"{i}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp", help="a tsv expression matrix (gene by sample) with row names and without col names.", type=str, required=True)
    parser.add_argument("--thr", help="get genes expressed in at least thr cells", type=int, required=True)
    parser.add_argument("--out", help="Path of the file which is a list of genes expressed", type=str, required=True)

    args = parser.parse_args()
    get_expressed_genes(args)