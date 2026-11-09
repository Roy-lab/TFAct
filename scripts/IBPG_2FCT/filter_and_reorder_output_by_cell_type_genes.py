import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("--indir", help="ibpg_2fct output directory", type=str, required=True)
parser.add_argument("--outdir", help="output directory for outputing reordered cell type specific matrices", type=str, required=True)
parser.add_argument("--in_gene_order", help="gene used by ibpg_2fct. It must be a superset of cell_type_gene_order", type=str, required=True)
parser.add_argument("--cell_type", help="cell cluster name for cell type specific matrices to be extracted", type=str, required=True)
parser.add_argument("--cell_type_gene_order", help="genes used for benchmarking", type=str, required=True)
parser.add_argument("--num_rep", help="number of stability rep", type=int, required=True)

print("Note: currently only 7 digits are saved for each entry to save space")

args = parser.parse_args()
indir = args.indir
outdir = args.outdir
in_gene_order = np.loadtxt(args.in_gene_order, dtype=str) 
cell_type = args.cell_type
cell_type_gene_order = np.loadtxt(args.cell_type_gene_order, dtype=str) 
num_rep = args.num_rep

gene_dict={in_gene_order[i]:i for i in range(len(in_gene_order))}

order = [gene_dict[g] for g in cell_type_gene_order]

for i in range(num_rep):
    A_in = indir + "/rep_" + str(i) + "/" + cell_type + "/A.txt"
    P_in = indir + "/rep_" + str(i) + "/" + cell_type + "/P.txt"
    A_out = outdir + "/rep_" + str(i) + "/A.txt"
    P_out = outdir + "/rep_" + str(i) + "/P.txt"

    os.system("mkdir -p "+ outdir + "/rep_" + str(i))

    A_in_arr = np.loadtxt(A_in)
    P_in_arr = np.loadtxt(P_in)
    relative_err = indir + "/rep_" + str(i) + "/relative_error.txt"
    out_relative_err = outdir + "/rep_" + str(i) + "/relative_error.txt"

    np.savetxt(P_out, P_in_arr, fmt='%.7f')
    np.savetxt(A_out, A_in_arr[order,], fmt='%.7f')
    os.system("cp "+ relative_err + " " + out_relative_err)
