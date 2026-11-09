
nseed=$1 #50
maxReg=$2 #50
cellfile=$3 #/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/celltype_order.txt
indir=$4 #/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/Results/scCVN/liger_4cells_sqrt_19genesrm_varthre0.05_k10/pg0.2_pm0.8_pr0.2_maxReg50_b0.9/subsample/
nogid=$5
deletefolder=$6
make_con_exe=${7:-/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/brand_multiome_project/scMTNI_with_TFA/scripts/Makeconsensusnetwork1.sh}
estEdge=${8:-/mnt/dv/wid/projects2/Roy-common/programs/programs/estimateedgeconf/estimateEdgeConf}
####################################################################################
# Step0: Merge separate networks from different gene sets
####################################################################################
for cell in `cat $cellfile`
do
    bash ${make_con_exe} $nseed $maxReg $cell $indir $nogid $estEdge &
done
wait


## Delete ogid sub directories
if [[ $deletefolder == 1 ]]
then
  bash /mnt/dv/wid/projects7/Roy-singlecell2/suvo_work/scMTNI_brandlab_human_Combined/inputs_for_GRN_without_TFA/scripts_b4_q1/deleteogidsubdirectory.sh $nseed $maxReg $cellfile $indir 
fi

