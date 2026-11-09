#!/bin/bash

nseed=$1 #50
maxReg=$2 #50
cell=$3 
indir=$4 #/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/brand_multiome_project/scMTNI_with_TFA/results/Lambda_0000_4_1
nogid=$5
####################################################################################
# Step0: Merge separate networks from different gene sets
####################################################################################
notdone=0
for r in `seq 0 1 $((nseed-1))` 
do
    {
    outdir=$indir/randseed${r}/$cell/fold0
    mkdir -p $outdir
    echo -n "" > $outdir/var_mb_pw_k${maxReg}.txt
    for og in `seq 0 1 $((nogid-1))`  
    do
        path0=$indir/randseed${r}/ogid${og}/
        path=${path0}/results/$cell/fold0
        #echo $path/var_mb_pw_k${maxReg}.txt 
        if [[ -s $path/var_mb_pw_k${maxReg}.txt ]];then
           cat $path/var_mb_pw_k${maxReg}.txt >> $outdir/var_mb_pw_k${maxReg}.txt
        else
           echo $path/var_mb_pw_k${maxReg}.txt "empty"
           notdone=1
        fi
    done
    #echo "Output: " $outdir
    #echo " "
    } &
done
wait

if [[ $notdone == 1 ]]
then 
    echo "Not all ogids are done"
    exit 0
fi
####################################################################################
# Step1: Make consensus network
####################################################################################

# Make a list of network files:
estEdge=${6:-/mnt/dv/wid/projects2/Roy-common/programs/programs/estimateedgeconf/estimateEdgeConf}
mkdir -p $indir/analysis/
cd $indir/analysis/
    echo $cell
    mkdir -p $cell
    echo -n "" > $cell/network_files.txt

    ## Prepare a list of all the subsample-specific network file paths for the given cell cluster
    ## The number of paths on the list will be equal to the number of subsamples.
    for rseed in `seq 0 1 $((nseed-1))`
    do
       ls ${indir}/randseed${rseed}/$cell/fold0/var_mb_pw_k${maxReg}.txt >> $cell/network_files.txt
    done

    ## Make a consensus network given the subsample-specific networks for the given cell cluster
    ${estEdge} $cell/network_files.txt 0 $cell/net_ alledges

## Remove the cell cluster suffix from the regulator and target gene names in the consensus network
sed "s/_${cell}//g;s/_${cell}//g" ${cell}/net_alledge.txt > ${cell}/consensus_edges.txt

