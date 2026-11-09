# indir should have ${species or cell cluster name}/nw.tsv and ${species or cell cluster name}/genes.txt
# will output to ${outdir}/${species or cell cluster name}_prior.txt and ${outdir}/lda_model_k${k}_${species or cell cluster name}.mat
# stdout would be logged at ${outdir}/${species or cell cluster name}.log
indir=$1 # /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/testoutput/0.2
outdir=$2 

script_dir=$(readlink -f ../../../scripts/)
prefixes=( $(ls ${indir}) )

tf_lst=${3:-/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt}
num_cells=${4:-2000}
cell_type=${5:-10} # cell state. We assume that each cell type has, say, 10 states, defined by cell type specific GRN and TFA profiles
num_rep=${6:-1}

for p in "${prefixes[@]}"; do
    {
    mkdir ${outdir}/${p}
    python ${script_dir}/generate_simulated_matrices_and_TFA_cell_topic.py --GRN ${indir}/${p}/nw.tsv \
	   --gene ${indir}/${p}/genes.txt \
	   --reg ${tf_lst} \
	   --num_cells ${num_cells} \
	   --cell_type ${cell_type} \
	   --num_rep ${num_rep} \
	   --outdir ${outdir}/${p}
    }&
done
wait
