#!/bin/bash

## The output dir must already exists

indir=$1 # contains output_randseed[0-9]+_[0-9]+.tar.gz
outdir=$2
num_subsamp=$3
num_og=$4
# logdir=$5 # contains randseed${subsamp}_${og}_scMTNI_Brand2024.out

## Create a dir for backing up the condor output tarballs
mkdir -p ${outdir}/tarballs

cpus=$(getconf _NPROCESSORS_ONLN)

function pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
        sleep 0.01
    done
}

## For each subsample
for subsamp in $(seq 0 $((${num_subsamp}-1)) );
do
   ## For each orthogroup
   {
   		for og in $(seq 0 $((${num_og}-1)) );
   		do 
			# This check nonexisting files. If this happens, change the queue file scMTNI_Brand2024_condor_queue_args.txt to
			# lines consisting of ${subsamp} ${og}, before condor_submit the corresponding submit file. We need to comment out
			# the code block after the "fi" and before the first "done" before running this bash script to check these files.
			if [ ! -f ${indir}/output_randseed${subsamp}_${og}.tar.gz ]; then
  				echo "${indir}/output_randseed${subsamp}_${og}.tar.gz does not exist"
			fi

			## Follow Shilu's preferred dir structure
			mkdir -p ${outdir}/randseed${subsamp}/ogid${og}/
			## Should extract a subdir named output_randseed${subsamp}_${og}/ 
			mv ${indir}/output_randseed${subsamp}_${og}.tar.gz ${outdir}/randseed${subsamp}/ogid${og}/results.tar.gz 
			tar -xzf ${outdir}/randseed${subsamp}/ogid${og}/results.tar.gz -C ${outdir}/randseed${subsamp}/ogid${og}/ --transform "s/output_randseed${subsamp}_${og}/results/"
   		done
   }&
   pwait $cpus
done
wait
echo "Finished"