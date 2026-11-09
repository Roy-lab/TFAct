
exp_tf_mat_dir=$1
cluster_file=$2
readarray -t cluster < $cluster_file
out_dir=${exp_tf_mat_dir}/filelists
mkdir -p ${out_dir}


outfile=${out_dir}/filelist.txt
cur_exp_dir=${exp_tf_mat_dir}
echo -n '' > ${outfile}
for c in "${cluster[@]}"; do
echo -e "${c}\t${cur_exp_dir}/${c}.table" >> ${outfile}
done

