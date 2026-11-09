#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data
# or at scMTNI_data/
# Also, need to specify a path to gsl lib
# export LD_LIBRARY_PATH=/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/brand_multiome_project/scMTNI_with_TFA/gsl/lib
#
outdir=TFA_output

ss_algos=( NetREX_CF HALS-p ibpg ibpg_2fct_2000_1000 ibpg_2fct_2000_100 ibpg_2fct_100_1000 ibpg_2fct kl_nmf_2fct kl_nmf_2fct_2000_1000 )

other_algos=( NINCA_quad TIGER regularized_nca supirfactor_hierarchical supirfactor_shallow )

algos=( NetREX_CF HALS-p ibpg ibpg_2fct_2000_1000 ibpg_2fct_2000_100 ibpg_2fct_100_1000 ibpg_2fct kl_nmf_2fct kl_nmf_2fct_2000_1000 NINCA_quad TIGER regularized_nca supirfactor_hierarchical supirfactor_shallow )

compare=/mnt/dv/wid/projects2/Roy-common/programs/programs/comparetargetsets_allnets_v3 #/mnt/dv/wid/projects2/Roy-regnet-inference/yeast_regnet/yeast_regnet_from_scratch/may_6/test_new_ptf/comparetargetsets_allnets/compareTgtSetAcrossDataset

mouse_gs=( mESC_chipunion_KDUnion_intersect mESC_KDUnion mESC_chipunion )

human_gs=( Cus_KO_0.01_Cus_ChIP_2_intersect.gs_network Cus_ChIP_2.gs_network Cus_KO_0.01.gs_network Cus_KO_0.01_union_intersect.gs_network UniBind.gs_network )

top_num=10000
top=_top10000
other_algo_top=_top_10000

# TFA_output/kl_nmf_2fct/processed_FBS_reprogramming/cluster1/A_top30000.

# TFA_output/regularized_nca/processed_Ranzoni_fetal_hematopoiesis_fine_tree/GPs_stability/Cus_ChIP_2.gs_network/stats_10000/Enr_filelist_top_10000.txt

# FBS
FBS_clusters=( cluster1 cluster2 cluster3 cluster4 cluster5 cluster7 cluster8 )
FBS_ds=processed_FBS_reprogramming
# BH
BH_clusters=( cluster1 cluster2 cluster3 cluster6 cluster7 cluster8 cluster9 cluster10 )
BH_ds=processed_Buenrostro_Hematopoiesis
# RFH
RFH_clusters=( B-cells Erythroid-cells GPs Granulocytes HSC-MPPs LMPs Mast-cells Megakaryocytes MEMPs Monocytes pDCs )
RFH_ds=processed_Ranzoni_fetal_hematopoiesis_fine_tree


# ss_algos
for c in "${FBS_clusters[@]}"; do
    {
    for a in "${ss_algos[@]}"; do
        for gs in "${mouse_gs[@]}"; do
            cur_path=${outdir}/${a}/${FBS_ds}/${c}/A${top}.${gs}.txt
            mkdir -p ${outdir}/${a}/${FBS_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/${a}/${FBS_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} 0.05 ${cur_out}
        done
    done
    } &
done
wait

# other algorithms
for c in "${FBS_clusters[@]}"; do
    {
    for a in "${other_algos[@]}"; do
        for gs in "${mouse_gs[@]}"; do
            cur_path=${outdir}/${a}/${FBS_ds}/${c}_stability/${gs}/stats_${top_num}/Enr_filelist${other_algo_top}.txt
            mkdir -p ${outdir}/${a}/${FBS_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/${a}/${FBS_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} 0.05 ${cur_out}
        done
    done
    } &
done
wait





# ss_algos
for c in "${BH_clusters[@]}"; do
    {
    for a in "${ss_algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/${a}/${BH_ds}/${c}/A${top}.${gs}.txt
            mkdir -p ${outdir}/${a}/${BH_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/${a}/${BH_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} 0.05 ${cur_out}
        done
    done
    } &
done
wait

# other algorithms
for c in "${BH_clusters[@]}"; do
    {
    for a in "${other_algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/${a}/${BH_ds}/${c}_stability/${gs}/stats_${top_num}/Enr_filelist${other_algo_top}.txt
            mkdir -p ${outdir}/${a}/${BH_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/${a}/${BH_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} 0.05 ${cur_out}
        done
    done
    } &
done
wait



# ss_algos
for c in "${RFH_clusters[@]}"; do
    {
    for a in "${ss_algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/${a}/${RFH_ds}/${c}/A${top}.${gs}.txt
            mkdir -p ${outdir}/${a}/${RFH_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/${a}/${RFH_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} 0.05 ${cur_out}
        done
    done
    } &
done
wait

# other algorithms
for c in "${RFH_clusters[@]}"; do
    {
    for a in "${other_algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/${a}/${RFH_ds}/${c}_stability/${gs}/stats_${top_num}/Enr_filelist${other_algo_top}.txt
            mkdir -p ${outdir}/${a}/${RFH_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/${a}/${RFH_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} 0.05 ${cur_out}
        done
    done
    } &
done
wait

## count
outfile_all=predictable_tf_nums.txt
>${outfile_all}

for c in "${FBS_clusters[@]}"; do
    for a in "${algos[@]}"; do
        for gs in "${mouse_gs[@]}"; do
            cur_path=${outdir}/${a}/${FBS_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "${a}\tFBS\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done
done

for c in "${BH_clusters[@]}"; do
    for a in "${algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/${a}/${BH_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "${a}\tBH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done
done


for c in "${RFH_clusters[@]}"; do
    for a in "${algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/${a}/${RFH_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "${a}\tRFH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done
done
