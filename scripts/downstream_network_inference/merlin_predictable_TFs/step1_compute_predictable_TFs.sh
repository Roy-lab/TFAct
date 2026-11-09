# TODO: add back other datasets
# need to specify a path to gsl lib
# export LD_LIBRARY_PATH=/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/brand_multiome_project/scMTNI_with_TFA/gsl/lib
#
outdir=GRN_outputs

algos=( no_TFA NetREX_CF ibpg_2fct ibpg_2fct_2000_1000 ibpg_2fct_2000_100 ibpg_2fct_100_1000 kl_nmf_2fct kl_nmf_2fct_2000_1000 NINCA_quad TIGER supirfactor_hierarchical supirfactor_shallow regularized_nca )

compare=/mnt/dv/wid/projects2/Roy-common/programs/programs/comparetargetsets_allnets_v3_nointersect/compareTgtSetAcrossDataset

mouse_gs=( mESC_chipunion_KDUnion_intersect mESC_KDUnion mESC_chipunion )

human_gs=( Cus_KO_0.01_Cus_ChIP_2_intersect.gs_network Cus_ChIP_2.gs_network Cus_KO_0.01.gs_network Cus_KO_0.01_union_intersect.gs_network UniBind.gs_network )

top_num=${1:-10000}
top=_top${top_num}
other_algo_top=_top_${top_num}

merlinc=5
merlinp=-5
merlinr=4
merlin_hyper=c${merlinc}_p${merlinp}_r${merlinr}

edgetype=( gene_only TFA_only max_conf )

gname=${outdir}/algos/goldnames.txt
echo "gt" > ${gname}

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

##
## Part 1: replaced
##
# 
# for c in "${FBS_clusters[@]}"; do
#     {
#     for a in "${algos[@]}"; do
#         for gs in "${mouse_gs[@]}"; do
#             cur_path=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${FBS_ds}/${c}/stability/${gs}/stats_${top_num}/Enr_filelist${other_algo_top}.txt
#             mkdir -p ${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${FBS_ds}/${c}/A${top}.${gs}
#             cur_out=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${FBS_ds}/${c}/A${top}.${gs}/
#             ${compare} ${cur_path} ${gname} 0.05 ${cur_out}
#         done
#     done
#     } &
# done
# wait


for c in "${BH_clusters[@]}"; do
    {
    for a in "${algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${BH_ds}/${c}/stability/${gs}/stats_${top_num}/Enr_filelist${other_algo_top}.txt
            mkdir -p ${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${BH_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${BH_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} ${gname} 0.05 ${cur_out}
        done
    done
    } &
done
wait


for c in "${RFH_clusters[@]}"; do
    {
    for a in "${algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${RFH_ds}/${c}/stability/${gs}/stats_${top_num}/Enr_filelist${other_algo_top}.txt
            mkdir -p ${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${RFH_ds}/${c}/A${top}.${gs}
            cur_out=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${RFH_ds}/${c}/A${top}.${gs}/
            ${compare} ${cur_path} ${gname} 0.05 ${cur_out}
        done
    done
    } &
done
wait

# count
mkdir -p merlin_predictable_tf/${merlin_hyper}
outfile_all=merlin_predictable_tf/${merlin_hyper}/merlin_predictable_tf_nums_top${top_num}.txt
> ${outfile_all}

for c in "${FBS_clusters[@]}"; do
    for a in "${algos[@]}"; do
        for gs in "${mouse_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${FBS_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "${a}\tFBS\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done
done

for c in "${BH_clusters[@]}"; do
    for a in "${algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${BH_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "${a}\tBH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done
done

for c in "${RFH_clusters[@]}"; do
    for a in "${algos[@]}"; do
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_${a}_replaced/${merlin_hyper}/${RFH_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "${a}\tRFH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done
done


###
# Part 1 ends
###

###
## Part 2: concat
###
mkdir -p merlin_predictable_tf/${merlin_hyper}/concat
for etype in "${edgetype[@]}"; do
    {
    # for c in "${FBS_clusters[@]}"; do
    #     {
    #     for a in "${algos[@]}"; do
    #         if [ "${a}" = "no_TFA" ]; then
    #             continue
    #         fi
    #         for gs in "${mouse_gs[@]}"; do
    #             cur_path=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${FBS_ds}/${c}/stability/${gs}/stats_${top_num}/${etype}_Enr_filelist${other_algo_top}.txt
    #             mkdir -p ${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${FBS_ds}/${c}/${etype}_A${top}.${gs}
    #             cur_out=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${FBS_ds}/${c}/${etype}_A${top}.${gs}/
    #             ${compare} ${cur_path} ${gname} 0.05 ${cur_out}
    #         done
    #     done
    #     } &
    # done
    # wait

    for c in "${BH_clusters[@]}"; do
        {
        for a in "${algos[@]}"; do
            if [ "${a}" = "no_TFA" ]; then
                continue
            fi
            for gs in "${human_gs[@]}"; do
                cur_path=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${BH_ds}/${c}/stability/${gs}/stats_${top_num}/${etype}_Enr_filelist${other_algo_top}.txt
                mkdir -p ${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${BH_ds}/${c}/${etype}_A${top}.${gs}
                cur_out=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${BH_ds}/${c}/${etype}_A${top}.${gs}/
                ${compare} ${cur_path} ${gname} 0.05 ${cur_out}
            done
        done
        } &
    done
    wait

    for c in "${RFH_clusters[@]}"; do
        {
        for a in "${algos[@]}"; do
            if [ "${a}" = "no_TFA" ]; then
                continue
            fi
            for gs in "${human_gs[@]}"; do
                cur_path=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${RFH_ds}/${c}/stability/${gs}/stats_${top_num}/${etype}_Enr_filelist${other_algo_top}.txt
                mkdir -p ${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${RFH_ds}/${c}/${etype}_A${top}.${gs}
                cur_out=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${RFH_ds}/${c}/${etype}_A${top}.${gs}/
                ${compare} ${cur_path} ${gname} 0.05 ${cur_out}
            done
        done
        } &
    done
    wait

    ## count
    outfile_all=merlin_predictable_tf/${merlin_hyper}/concat/${etype}_merlin_predictable_tf_nums_top${top_num}.txt
    > ${outfile_all}

    for c in "${FBS_clusters[@]}"; do
        for a in "${algos[@]}"; do
            if [ "${a}" = "no_TFA" ]; then
                continue
            fi
            for gs in "${mouse_gs[@]}"; do
                cur_path=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${FBS_ds}/${c}/${etype}_A${top}.${gs}/_factor_coenr.txt
                cnt=`cat $cur_path | wc -l`
                cnt=`expr $cnt - 1`
                echo -e "${a}\tFBS\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
            done
        done
        # add no_TFA
        for gs in "${mouse_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_no_TFA_replaced/${merlin_hyper}/${FBS_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "no_TFA\tFBS\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done

    for c in "${BH_clusters[@]}"; do
        for a in "${algos[@]}"; do
            if [ "${a}" = "no_TFA" ]; then
                continue
            fi
            for gs in "${human_gs[@]}"; do
                cur_path=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${BH_ds}/${c}/${etype}_A${top}.${gs}/_factor_coenr.txt
                cnt=`cat $cur_path | wc -l`
                cnt=`expr $cnt - 1`
                echo -e "${a}\tBH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
            done
        done
        # add no_TFA
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_no_TFA_replaced/${merlin_hyper}/${BH_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "no_TFA\tBH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done

    for c in "${RFH_clusters[@]}"; do
        for a in "${algos[@]}"; do
            if [ "${a}" = "no_TFA" ]; then
                continue
            fi
            for gs in "${human_gs[@]}"; do
                cur_path=${outdir}/algos/merlin_${a}_concat/${merlin_hyper}/${RFH_ds}/${c}/${etype}_A${top}.${gs}/_factor_coenr.txt
                cnt=`cat $cur_path | wc -l`
                cnt=`expr $cnt - 1`
                echo -e "${a}\tRFH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
            done
        done
        # add no_TFA
        for gs in "${human_gs[@]}"; do
            cur_path=${outdir}/algos/merlin_no_TFA_replaced/${merlin_hyper}/${RFH_ds}/${c}/A${top}.${gs}/_factor_coenr.txt
            cnt=`cat $cur_path | wc -l`
            cnt=`expr $cnt - 1`
            echo -e "no_TFA\tRFH\t${c}\t${gs}\t${cnt}" >> ${outfile_all}
        done
    done

    } &
done
wait
