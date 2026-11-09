# indir should have ${species or cell cluster name}/nw.tsv and ${species or cell cluster name}/genes.txt
# will output to ${outdir}/${species or cell cluster name}_prior.txt and ${outdir}/lda_model_k${k}_${species or cell cluster name}.mat
# stdout would be logged at ${outdir}/${species or cell cluster name}.log
indir=$1 # /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/testoutput/0.2
outdir=$2 
k=$3 # number of topics
multiplicity=$4

script_dir=$(readlink -f ../../../scripts/)
prefixes=( $(ls ${indir}) )

tf_lst=${5:-/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt}

for p in "${prefixes[@]}"; do
    matlab -nodisplay -nodesktop -r "addpath('${script_dir}'); add_noise_lda('${indir}/${p}/nw.tsv', ${k}, ${multiplicity}, '${p}', '${outdir}'); quit" \
       >  ${outdir}/${p}.log &
done
wait

