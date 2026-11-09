# assume that we have already put the cell lineage tree in outdir (though it is not required for running the script)

in_gene_dir=$1 # /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/testoutput/0.2
in_exp_prior_dir=$2 # /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/inputtoTFA/0.2/k10
outdir=$3 # /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/inputtoTFA/0.2/k10/all_exps

#TFA_sim_rep=( rep0 )
#exp_name=( E E_0.3 E_0.6 )
#cell_types=( Anc2 Anc5 Calb Cgla Sbay Spom right )

IFS=',' read -r -a TFA_sim_rep <<< $4
IFS=',' read -r -a exp_name <<< $5
IFS=',' read -r -a cell_types <<< $6
IFS=',' read -r -a cell_types_w_prior_available <<< $7

# the below two should correspond to each other, and cell_types_w_prior_available should be a subset of cell_types
#cell_types_w_prior_available=( Anc2 Anc5 Calb Cgla Sbay Spom right )
cell_types_prior_confidence=( 1 1 1 1 1 1 1 )

##
## the file storing paths to expression files will be put as ${outdir}/${rep}/${exp_n}_all_exps.txt
## the file storing paths to priors along with their confidence will be put as ${outdir}/all_priors.txt
## the file storing the union of genes is stored as ${outdir}/all_genes.txt
##

## make union of genes
mkdir -p ${outdir}
gene_files=()
for c in "${!cell_types[@]}"; do
    gene_files[c]=${in_gene_dir}/${cell_types[c]}/genes.txt
done

cat ${gene_files[@]} | sort -u > ${outdir}/all_genes.txt

# now make the file containing expression matrices' paths
for rep in "${TFA_sim_rep[@]}"; do
    for E in "${exp_name[@]}"; do
        mkdir -p ${outdir}/${rep}
        out_exp_file=${outdir}/${rep}/${E}_all_exps.txt
        echo -n "" > ${out_exp_file}
        for c in "${cell_types[@]}"; do
            echo -e "${c}\t${in_exp_prior_dir}/${c}/${rep}/${E}.txt" >> ${out_exp_file}
        done
    done
done

# now make the file containing priors' paths
out_prior_file=${outdir}/all_priors.txt
echo -n "" > ${out_prior_file}
for i in "${!cell_types[@]}"; do
    echo -e "${cell_types[i]}\t${in_exp_prior_dir}/${cell_types[i]}_prior.txt\t${cell_types_prior_confidence[i]}" >> ${out_prior_file}
done
