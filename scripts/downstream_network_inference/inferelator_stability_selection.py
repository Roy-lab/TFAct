import argparse
import numpy as np
import gzip
from collections import defaultdict

# GRN_outputs/algos/inferelator_default/processed_FBS_reprogramming/cluster1/4/network.tsv.gz
def stability_selection(args):
    edge_dict = defaultdict(int)
    for i in range(args.num_subsamp):
        cur_file = args.indir + '/' + str(i) + '/network.tsv.gz'
        with gzip.open(cur_file, 'rt') as f:
            next(f)
            for line in f:
                l = line.strip().split('\t')
                edge_dict[(l[1],l[0])]+=1 # target\tregulator
    for k in edge_dict:
        edge_dict[k] /= args.num_subsamp
    edge_lst = sorted(edge_dict.items(), key=lambda x: x[1], reverse=True)
    with open(args.outfile, 'w') as f:
        for k, v in edge_lst:
            dummy = f.write(k[0] + '\t' + k[1] + '\t' + str(v)+'\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", help="directory that contains (subsample)/network.tsv.gz. The file's first two fields are target\\tregulator.", type=str, required=True)
    parser.add_argument("--num_subsamp", help='Number of subsamples', type=int, required=True)
    parser.add_argument("--outfile", help='output tsv file name. It will output regulator\\ttarget\\tconfidence', type=str, required=True)

    args = parser.parse_args()
    
    stability_selection(args)