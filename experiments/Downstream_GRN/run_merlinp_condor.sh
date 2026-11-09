#!/bin/sh

# ## Exit immediately if a command exits with a non-zero status
# set -e

# "$(subsamp) $(clust) 5 -5 4"

subsamp=$1
clust=$2
prior_conf=$3
p=$4
r=$5
dataset=$6
algo=$7
ds_action=$8

tar -xzf libgsl.tar.gz --transform "s/lib/libgsl/"
export LD_LIBRARY_PATH=$(pwd)/libgsl

echo -e "${dataset}_${clust}_${subsamp}\t${clust}_network_TFA.txt\t${prior_conf}" > condor_config.txt

outdir=output_${dataset}_${subsamp}_${clust}_${algo}_TFA${ds_action}

if [ -f $outdir.tar.gz ]
then
    tar xzf $outdir.tar.gz
else
    mkdir $outdir
    tar czf $outdir.tar.gz $outdir
fi

clean_up() {
  echo "Job is exiting; tar result"
  tar -czf $outdir.tar.gz $outdir
  rm -rf $outdir
  rm -rf libgsl*
  rm condor_config.txt
  exit $1
}

trap clean_up SIGTERM SIGABRT

if [ -f $outdir/fold0/prediction_k300.txt ]
then
    echo "./merlin -d ${subsamp}.txt -l ${clust}_complete_allregulators_TFA.txt -o ${outdir} -q condor_config.txt -p $p -r $r -e 1 -a 1"
    ./merlin -d ${subsamp}.txt -l ${clust}_complete_allregulators_TFA.txt -o ${outdir} -q condor_config.txt -p $p -r $r -e 1 -a 1 &
    PID=$!
    wait $PID
    EXIT_STATUS=$?
else
    echo "./merlin -d ${subsamp}.txt -l ${clust}_complete_allregulators_TFA.txt -o ${outdir} -q condor_config.txt -p $p -r $r -e 1"
    ./merlin -d ${subsamp}.txt -l ${clust}_complete_allregulators_TFA.txt -o ${outdir} -q condor_config.txt -p $p -r $r -e 1 &
    PID=$!
    wait $PID
    EXIT_STATUS=$?
fi


# cleanup                                                                                                                                                                                                                                                                                                                                    
# check the exit code of our merlin execute                                                                                                                                                                                                                                                                                                  
if [ $EXIT_STATUS -ne 0 ]
then
  echo "there was an error: $EXIT_STATUS"
fi

clean_up $EXIT_STATUS



