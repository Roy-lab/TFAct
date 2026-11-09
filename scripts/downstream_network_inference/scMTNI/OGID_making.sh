COUNTER=1
cluster_file=$1 # /mnt/dv/wid/projects7/Roy-singlecell2/suvo_work/leukemic_patient_specific_study/Patient1_BM4/scMTNI/celltype_order.txt
readarray -t cluster < $cluster_file

echo -n "" > OGID.txt

cat Union.txt | while read line

do

cur_line="OG${COUNTER}_0\t"
for c in "${cluster[@]}"; do
    cur_line+="${line}_${c},"
done

cur_line="${cur_line%,}"

echo -e "${cur_line}" >> OGID.txt

echo $COUNTER

COUNTER=$((COUNTER+1))

done

 

sed -i '1 i\orthoID\tgenes' OGID.txt
