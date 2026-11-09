import os
import pandas as pd
import argparse

#filelist='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/gammlab/data/simplenmf_union_randnmf_ligernorm_k20_iter5/filelistbycluster.txt'
#indir='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/gammlab/data/simplenmf_union_randnmf_ligernorm_k20_iter5/'
def main(args):
    indir=args.indir+'/'
    nseed=args.nseed
    filelist=pd.read_table(args.filelist,header=None)
    infiles=filelist[1].values
    cells=filelist[0].values
    frac=2/3 ##args.fraction
    for seed in range(nseed):
        merged=[]
        for i in range(len(cells)):
            cell=cells[i]
            print(infiles[i])
            data = pd.read_table(infiles[i], sep='\t',index_col=0)
            data1=data.transpose()
            data2=data1.sample(n=int(data1.shape[0]*frac), replace=False,random_state=seed)  ##12067
            data3=data2.transpose()
            outdir=indir+'/subsample/randseed'+str(seed)+'/'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            outname=outdir+cell+'.table'
            data3.to_csv(outname, index = True, index_label='Gene', header=True, sep='\t', mode='w',float_format='%g')
            data3.index=data3.index.str.replace('_'+cell,'')
            merged.append(data3)
            print(data3.shape)

        #dfmerge = pd.concat(merged, axis=1)
        #outname1=outdir+'merged.table'
        #dfmerge.to_csv(outname1, index = True, index_label='Gene', header=True, sep='\t', mode='w',float_format='%g')
        #print(dfmerge.shape)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--filelist',

                        help='filelist',

                        type=str,

                        default='')
    parser.add_argument('--nseed',

                        help='number of subsamples',

                        type=int,

                        default=50)   
    parser.add_argument('--indir',

                        help='filelist',

                        type=str,

                        default='')
    parser.add_argument('--fraction',

                        help='fraction of samples',

                        type=float,

                        default=0.5)

    args = parser.parse_args()
    main(args)



