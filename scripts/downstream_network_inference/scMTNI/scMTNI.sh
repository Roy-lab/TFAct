
indir=$1 
outdir=$2

filelist_dir=$3 
regfile=$4
exe=$5


cur_outdir=${outdir}/results
mkdir -p ${cur_outdir}
filelist=${filelist_dir}/filelist.txt
python ${exe} --filelist $filelist --regfile $regfile --indir ${indir}/ --outdir ${cur_outdir}/ --splitgene 50 --motifs 1

