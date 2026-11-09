import os
import numpy as np
import argparse
import glob

parser = argparse.ArgumentParser()
# We assume that each directory containing A.txt and P.txt also has a relative_error.txt, where the last line indicates the final relative error.
parser.add_argument("--A_file_pattern", help="the regex that matches all files for matrix A (gene by tf, all numeric and tab separated). Please enclose with single quotes.", type=str, required=True)
parser.add_argument("--P_file_pattern", help="the regex that matches all files for matrix P (tf by cell, all numeric and tab separated). Please enclose with single quotes.", type=str, required=True)

parser.add_argument("--frac_ignored", help="From 0 to 1. ignore this fraction of matrices corresponding to top relative loss", type=float, default=0)

parser.add_argument("--output_precision", help="The precision of entries of output files", type=int, default=8)
parser.add_argument("--out_A", help="The output matrix for A averaged over GRN matrices that pass the filter.", type=str, required=True)
parser.add_argument("--out_P", help="The output matrix for P averaged over GRN matrices that pass the filter", type=str, required=True)

args = parser.parse_args()

file_pattern = [args.A_file_pattern, args.P_file_pattern]
outfile = [args.out_A, args.out_P]
for i in range(len(file_pattern)):
    file_list = glob.glob(file_pattern[i])
    rel_file_list = [os.path.dirname(_)+ '/relative_error.txt' for _ in file_list]
    rel_err = []
    for rel_file in rel_file_list:
        with open(rel_file, 'rb') as f:
            try:  # catch OSError in case of a one line file 
                f.seek(-2, os.SEEK_END)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
            except OSError:
                f.seek(0)
            rel_err += [float(f.readline().decode().rstrip("\n "))]
    indices_rm=sorted(range(len(rel_err)), key=lambda i: rel_err[i])[-int(len(rel_err)*args.frac_ignored):]
    file_list = [file_list[_] for _ in range(len(file_list)) if _ not in indices_rm]
    matrix_sum = None
    file_count = 0
    for file_path in file_list:
        current_matrix = np.loadtxt(file_path)
        if matrix_sum is None:
            matrix_sum = current_matrix
        else:
            matrix_sum += current_matrix
        file_count += 1
    average_matrix = matrix_sum / file_count
    np.savetxt(outfile[i], average_matrix, fmt=f'%.{args.output_precision}f')
