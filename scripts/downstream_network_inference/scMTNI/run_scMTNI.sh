#!/bin/sh

## Exit immediately if a command exits with a non-zero status
set -e

## Copied from Shilu's script:
## /mnt/dv/wid/projects2/Roy-regnet-inference/CVN/CVNvrConstraint_SingleCellwSparsitywMotif_upto1transition/runCVN_condor.sh
##

subsamp=$1
og=$2

ogid=$3
proot=$4
inputTarget=$5  ## AllGenes.txt
maxReg=$6
b1=$7
b2=$8
tar -xzf prior.tar.gz

echo "Processing ${subsamp} og${og}..."
 
## Create cell-cluster specific output dirs
for folder in `cut -f 3 condor_config.txt`
do
    mkdir -p $folder
done


tar -xzf libgsl.tar.gz --transform "s/lib/libgsl/"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)/libgsl

## add -c yes for single cell data to split the AllGenes.txt
echo ./scMTNI -f condor_config.txt -x${maxReg} -v1 -l TFs_OGs.txt -n $inputTarget -d celltype_tree_ancestor.txt -m $ogid -s celltype_order.txt -p $proot -c yes -b $b1 -q $b2
./CVN -f condor_config.txt -x${maxReg} -v1 -l TFs_OGs.txt -n $inputTarget -d celltype_tree_ancestor.txt -m $ogid -s celltype_order.txt -p $proot -c yes -b $b1 -q $b2
# ./scMTNI -f config_condor.txt -x50 -v1 -l TFs_OGs.txt -n AllGenes0.txt -d celltype_tree_ancestor.txt -m testdata_ogids.txt -s celltype_order.txt -p 0.2 -c yes -b -3 -q 3
# ./scMTNI -f config_condor.txt -x${maxReg} -v1 -l TFs_OGs.txt -n $inputTarget -d celltype_tree_ancestor.txt -m $ogid -s celltype_order.txt -p $proot -c yes -b $b1 -q $b2

## The desired name of the output dir
outdir=output_${subsamp}_${og}

## Rename the output dir with the desired name
mv results ${outdir}

## Compress the output dir.
## To unzip, use tar -xzvf <>.tar.gz
tar -czvf ${outdir}.tar.gz ${outdir}

