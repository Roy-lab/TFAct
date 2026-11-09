import argparse
import itertools
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('outpath')
    parser.add_argument('outreg')
    parser.add_argument('outgene')
    parser.add_argument('--num-reg', type = int, default = 300)
    parser.add_argument('--num-genes', type = int, default = 1000, help='must be not smaller than --num-reg')
    parser.add_argument('--density', type = float, default = 0.4)
    args = parser.parse_args()

    n_reg = args.num_reg
    n_genes = args.num_genes
    density = args.density

    regulators = ['Y'+ str(i) for i in range(n_reg)]
    targets = ['Y'+ str(i) for i in range(n_genes)]
    ind = np.random.choice(n_genes * n_reg, size=int(density * n_genes * n_reg), replace=False).tolist()
    all_edges = [_ for _ in itertools.product(regulators, targets)]
    edges_chosen = [all_edges[_] for _ in ind]
    with open(args.outpath, 'w') as f:
        for _ in edges_chosen:
            dummy = f.write(_[0] + '\t' + _[1] + '\n')
    np.savetxt(args.outreg, np.unique([_[0] for _ in edges_chosen]), fmt='%s')
    np.savetxt(args.outgene, np.unique([_[1] for _ in edges_chosen]), fmt='%s')
    