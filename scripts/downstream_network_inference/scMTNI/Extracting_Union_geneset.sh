outdir=$1
cluster_file=$2

readarray -t cluster < $cluster_file

for cell in "${cluster[@]}";
do
    {
        awk 'NR>1 {print $1}' ${outdir}/${cell}.table | sed 's/_'"$cell"'$//g' > ${outdir}/${cell}_genelist.txt
    } &
done
wait

cat ${outdir}/*_genelist.txt | sort | uniq > ${outdir}/Union.txt
rm ${outdir}/*_genelist.txt

wait