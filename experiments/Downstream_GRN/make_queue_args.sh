# as we may not submit that many jobs at one time, we split by each dataset. We have four blocks.

run_num_subsamps=50
FBS_clusters=(cluster1 cluster2 cluster3 cluster4 cluster5 cluster7 cluster8)
BH_clusters=(cluster1 cluster2 cluster3 cluster6 cluster7 cluster8  cluster9  cluster10)
RFH_fine_clusters=(B-cells Erythroid-cells GPs Granulocytes HSC-MPPs LMPs Mast-cells Megakaryocytes MEMPs Monocytes pDCs)

algo=(NINCA_quad ibpg_2fct_100_1000 kl_nmf_2fct_2000_1000 TIGER ibpg_2fct_2000_100 ibpg_2fct_2000_1000 supirfactor_hierarchical supirfactor_shallow ibpg_2fct kl_nmf_2fct regularized_nca NetREX_CF)

ds_action=(concat replaced)
datasets=(processed_FBS_reprogramming processed_Buenrostro_Hematopoiesis processed_Ranzoni_fetal_hematopoiesis_fine_tree)
outfile=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/queue_args.txt


####
## Below, only one block can be uncommented at the same time
####

### block 1 starts ###

# > ${outfile}

# for alg in "${algo[@]}"; do
#     for act in "${ds_action[@]}"; do
#         for r in $(seq 0 $(("${run_num_subsamps}" - 1))); do
#             for c in "${FBS_clusters[@]}"; do
#                 echo "processed_FBS_reprogramming ${r} ${c} ${alg} ${act}" >> ${outfile}
#             done
#         done
#     done
# done

### block 1 ends

### block 2 starts ###

# > ${outfile}
# for alg in "${algo[@]}"; do
#     for act in "${ds_action[@]}"; do
#         for r in $(seq 0 $(("${run_num_subsamps}" - 1))); do
#             for c in "${BH_clusters[@]}"; do
#                 echo "processed_Buenrostro_Hematopoiesis ${r} ${c} ${alg} ${act}" >> ${outfile}
#             done
#         done
#     done
# done

# ### block 2 ends

# ### block 3 starts

# > ${outfile}
# for alg in "${algo[@]}"; do
#     for act in "${ds_action[@]}"; do
#         for r in $(seq 0 $(("${run_num_subsamps}" - 1))); do
#             for c in "${RFH_fine_clusters[@]}"; do
#                 echo "processed_Ranzoni_fetal_hematopoiesis_fine_tree ${r} ${c} ${alg} ${act}" >> ${outfile}
#             done
#         done
#     done
# done

# ### block 3 ends

# ### block 4 starts

# > ${outfile}
# for alg in "no_TFA"; do
#     for act in "replaced"; do
#         for r in $(seq 0 $(("${run_num_subsamps}" - 1))); do
#             for c in "${FBS_clusters[@]}"; do
#                 echo "processed_FBS_reprogramming ${r} ${c} ${alg} ${act}" >> ${outfile}
#             done
#             for c in "${BH_clusters[@]}"; do
#                 echo "processed_Buenrostro_Hematopoiesis ${r} ${c} ${alg} ${act}" >> ${outfile}
#             done
#             for c in "${RFH_fine_clusters[@]}"; do
#                 echo "processed_Ranzoni_fetal_hematopoiesis_fine_tree ${r} ${c} ${alg} ${act}" >> ${outfile}
#             done
#         done
#     done
# done

# ### block 4 ends