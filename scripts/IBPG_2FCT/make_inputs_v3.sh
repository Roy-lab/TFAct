# assume that we have already put the cell lineage tree in outdir (though it is not required for running the script)

in_dir=$1 
outdir=$2 

readarray -t cell_types < $3


#cell_types_w_prior_available=
cell_types_prior_confidence=1

##
## the file storing paths to expression files will be put as ${outdir}/${rep}/${exp_n}_all_exps.txt
## the file storing paths to priors along with their confidence will be put as ${outdir}/all_priors.txt
## the file storing the union of genes is stored as ${outdir}/all_genes.txt
##

## make union of genes
mkdir -p ${outdir}
gene_files=()
for c in "${!cell_types[@]}"; do
    echo ${cell_types[c]}
    gene_files[c]=${in_dir}/${cell_types[c]}_allGenes.txt
done

cat ${gene_files[@]} | sort -u > ${outdir}/all_genes.txt

# now make the file containing expression matrices' paths

out_exp_file=${outdir}/all_exps.txt
echo -n "" > ${out_exp_file}
for c in "${cell_types[@]}"; do
    echo -e "${c}\t${in_dir}/${c}.table" >> ${out_exp_file}
done

# now make the file containing priors' paths
out_prior_file=${outdir}/all_priors.txt
echo -n "" > ${out_prior_file}
for i in "${!cell_types[@]}"; do
    echo -e "${cell_types[i]}\t${in_dir}/${cell_types[i]}_network.txt\t${cell_types_prior_confidence}" >> ${out_prior_file}
done
