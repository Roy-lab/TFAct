### to organize this, one needs to 1. specify the indir (where the remapped .tar.gz condor output files can be found)
### and 2. provide an output directory. See "#change this" below for other specifications. It will give a queue arg for each 
### of the datasets of failed jobs, which can then be rerun using condor, and the log files are archived in ${failed_job_dir}/logs

run_num_subsamps=50
FBS_clusters=(cluster1 cluster2 cluster3 cluster4 cluster5 cluster7 cluster8)
BH_clusters=(cluster1 cluster2 cluster3 cluster6 cluster7 cluster8  cluster9  cluster10)
RFH_fine_clusters=(B-cells Erythroid-cells GPs Granulocytes HSC-MPPs LMPs Mast-cells Megakaryocytes MEMPs Monocytes pDCs)

algo=(NINCA_quad ibpg_2fct_100_1000 kl_nmf_2fct_2000_1000 TIGER ibpg_2fct_2000_100 ibpg_2fct_2000_1000 supirfactor_hierarchical supirfactor_shallow ibpg_2fct kl_nmf_2fct regularized_nca)

ds_action=(concat replaced)

# merlin hyperparams
cmerlin=5
pmerlin=-5
rmerlin=4

# change this
indir=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/results/merlinp_${cmerlin}_${pmerlin}_${rmerlin}
outdir=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/GRN_outputs/algos
logdir=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/logs/merlinp_5_-5_4 # where things are logged in the merlin submit file
failed_job_dir=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/GRN_outputs/failed_merlin_jobs # if it returns scrambled output network. There could be other errors
mkdir -p ${failed_job_dir}/logs

# FBS reprogramming
cur_dataset=processed_FBS_reprogramming

for alg in "${algo[@]}"; do
    for act in "${ds_action[@]}"; do
        {
        for r in $(seq 0 $(("${run_num_subsamps}" - 1))); do
            for c in "${FBS_clusters[@]}"; do
                if [ ! -f ${indir}/output_${cur_dataset}_${r}_${c}_${alg}_TFA${act}.tar.gz ]; then
                    echo "${indir}/output_${cur_dataset}_${r}_${c}_${alg}_TFA${act}.tar.gz does not exist"
                    continue
                fi
                ### echo "${indir}/output_${cur_dataset}_${r}_${c}_${alg}_TFA${act}.tar.gz"
                cur_clust_outdir=${outdir}/merlin_${alg}_${act}/c${cmerlin}_p${pmerlin}_r${rmerlin}/${cur_dataset}/${c}
                mkdir -p ${cur_clust_outdir}
                ## Should extract a subdir named output_randseed${subsamp}_${og}/                                                                                                                             
                mv ${indir}/output_${cur_dataset}_${r}_${c}_${alg}_TFA${act}.tar.gz ${cur_clust_outdir}/${r}.tar.gz
                tar -xzf ${cur_clust_outdir}/${r}.tar.gz -C ${cur_clust_outdir}/ --transform "s/output_${cur_dataset}_${r}_${c}_${alg}_TFA${act}/${r}/"

                # check if the network has the wrong format (chtc issue)
                awk -F'\t' 'NF!=3{printf "\n"; exit}' ${cur_clust_outdir}/${r}/fold0/prediction_k300.txt > ${cur_clust_outdir}/${r}/fold0/k300_scrambled.txt.err
            done
        done
        }&
    done
done
wait

## we need to ensure that the following two files are empty after running the below. If not, copy the queue args to the one queued in the 
## merlin submit script (i.e. /mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/submit_files/submit_merlinp_p-5_r4.sub) and submit to condor
> ${failed_job_dir}/${cur_dataset}_scrambled_networks_logfiles.txt
> ${failed_job_dir}/${cur_dataset}_queue_args.txt
# compile failed jobs
for alg in "${algo[@]}"; do
    for act in "${ds_action[@]}"; do
        for r in $(seq 0 $(("${run_num_subsamps}" - 1))); do
            for c in "${FBS_clusters[@]}"; do
                cur_clust_outdir=${outdir}/merlin_${alg}_${act}/c${cmerlin}_p${pmerlin}_r${rmerlin}/${cur_dataset}/${c}
                if [ "$(wc -l < ${cur_clust_outdir}/${r}/fold0/k300_scrambled.txt.err)" -gt 0 ]; then
                    echo "${logdir}/${cur_dataset}_${alg}_${act}_${r}_${c}.log" >> ${failed_job_dir}/${cur_dataset}_scrambled_networks_logfiles.txt
                    echo "${cur_dataset} ${r} ${c} ${alg} ${act}" >> ${failed_job_dir}/${cur_dataset}_queue_args.txt
                    cp ${logdir}/${cur_dataset}_${alg}_${act}_${r}_${c}.* ${failed_job_dir}/logs
                fi
            done
        done
    done
done
wait

