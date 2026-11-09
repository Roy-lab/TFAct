import argparse
import numpy as np
import pandas as pd

def convert_prior(args): 
    gene=np.loadtxt(args.gene, dtype=str) # "testoutput/0.2/Anc5/genes.txt"
    tfs=np.loadtxt(args.tf_lst, dtype=str)
    prior_A_lst = np.loadtxt(args.prior_A, dtype=str) 
    gene_dict={gene[i]:i for i in range(len(gene))}
    reg_dict={tfs[i]:i for i in range(len(tfs))}

    num_gene=len(gene)
    num_tfs = len(tfs)

    priorA = np.zeros((num_gene,num_tfs))

    for r, t, _ in prior_A_lst:
        if r in tfs and t in gene:
            priorA[gene_dict[t], reg_dict[r]]=_

    df = pd.DataFrame(data=priorA,    
              index=gene,    
              columns=tfs)
    df.to_csv(args.out, sep='\t', index=True, compression='gzip')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene", help="a column file where each line is the name of a gene. This specifies the order of rows.", type=str, required=True)
    parser.add_argument("--tf_lst", help="a column file where each line is the name of a regulator. This specifies the order of cols.", type=str, required=True)
    parser.add_argument("--prior_A", help="The prior GRN edge file used where each line is regulator[tab]target[tab]confidence.", type=str, required=True)
    parser.add_argument("--out", help="Path of the compressed prior file, as required by inferelator.", type=str, required=True)
    
    args = parser.parse_args()
    convert_prior(args)