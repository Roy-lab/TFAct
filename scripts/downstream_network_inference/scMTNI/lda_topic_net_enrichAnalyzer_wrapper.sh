###
## Note that LDA is 1-based, in contrast to DTM, which is 0-based.
# The step 4: combine heatmaps below are commented out. The heatmap can be huge and not readily interpretable.
##
###
## This script has 3 configuration options: if all_regs_included is T, then all regulators (including _TFA) are considered for enrichment analysis, along with all targets. 
# If target_only is T given that all_regs_included is F, then only the target genes are considered for enrichment analysis. If both are F (default), then both mRNA-based regulators 
# (excluding _TFA regulators by downstream processing) and targets are considered for enrichment analysis.
###
result_dir=$1 # /mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/dynamic_networks/Andrew_multiome_project/full_tree_out/k10
lda_k=$2
go_regnet=${3:-/mnt/dv/wid/projects2/Roy-common/data/data_new/mouse/go/mm10/mm10_goterms_regnet.txt}
bg_genelist=${4:-/mnt/dv/wid/projects6/Roy-singlecell3/sridharanlab2/Andrew_multiome/spencer_work/results/all_genes.txt}
cell_type_order=${5:-/mnt/dv/wid/projects6/Roy-singlecell3/sridharanlab2/Andrew_multiome/spencer_work/results/scMTNI_in/celltype_order_viz.txt}
target_only=${6:-F}
indir=${7:-${result_dir}}
all_regs_included=${8:-F} # also consider _TFA regulators. This overwrite the target_only option.
enrichAnalyzerscript=${9:-/mnt/dv/wid/projects2/Roy-common/programs/programs/enrichanalyzer_Nongraph_Qval/enrichAnalyzer}

###
## Step 1: make input files for enrichAnalyzer
###
for i in $(seq 1 ${lda_k}); do 
    lda_result_dir=${result_dir}
    mkdir -p ${lda_result_dir}/ldaGOenrichment
    
    if [[ "$all_regs_included" == "T" ]]; then
        tail -n +2 ${indir}/networks_lda_k${lda_k}_topic${i}.txt \
            | awk -F'\t' '{
                reg=$1
                sub(/_TFA$/, "", reg)
                print reg"\t"$4
                print $2"\t"$4
            }' \
            | sort -u > ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset.txt

    elif [[ "$target_only" == "T" ]]; then
        tail -n +2 ${indir}/networks_lda_k${lda_k}_topic${i}.txt \
            | awk -F'\t' '{print $2"\t"$4}' \
            | sort -u > ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset.txt
    else
        tail -n +2 ${indir}/networks_lda_k${lda_k}_topic${i}.txt \
            | awk -F'\t' '{print $1"\t"$4; print $2"\t"$4}' \
            | sort -u > ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset.txt
    fi
    
    python2.7 /mnt/dv/wid/projects7/Roy-singlecell2/sahalab/scripts/makeGroup.py \
        ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset.txt \
        ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset_enrichAnalyzer.txt 5
    
    sed 's/^Cluster\([^[:space:]]\+\)/\1/' ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset_enrichAnalyzer.txt \
        > ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset_enrichAnalyzer.txt.tmp \
        && mv ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset_enrichAnalyzer.txt.tmp ${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset_enrichAnalyzer.txt
done

###
## Step 2: run enrichment
###
for i in $(seq 1 ${lda_k}); do 
{
    lda_result_dir=${result_dir}
    enrichAnalyzerInput=${lda_result_dir}/ldaGOenrichment/networks_lda_k${lda_k}_topic${i}_geneset_enrichAnalyzer.txt
    outdir=${lda_result_dir}/ldaGOenrichment/GOenrichment
    mkdir -p $outdir
    ${enrichAnalyzerscript} "${enrichAnalyzerInput}" "${bg_genelist}" "${go_regnet}" 0.05 "${outdir}/go_enrichAnalysis_topic${i}" persg
} &
done
wait

###
## Step 3: make heatmaps
###

script=/mnt/dv/wid/projects5/Roy-singlecell/js_work/scripts/enrichment_test_heatmapawk/enrichment_details_all_terms_to_hmawk.pl
hm_script=/mnt/dv/wid/projects5/Roy-singlecell/js_work/scripts/enrichment_test_heatmapawk/Heatmap.awk
W=20            # width of heatmap cell                                                                                                                                               
H=12            # height of heatmap cell 
MAX=10          # maximum of enrichment (-logP)

for i in $(seq 1 ${lda_k}); do
    lda_result_dir=${result_dir}
    GO_dir=${lda_result_dir}/ldaGOenrichment/GOenrichment
    DETAILS=${GO_dir}/go_enrichAnalysis_topic${i}_details.txt
    LIST=${cell_type_order}
    HMAWK="${DETAILS%.*}".hmawk

    ${script} ${DETAILS} ${LIST} $H $W > $HMAWK

    OUTFN="${HMAWK%.*}"  
    # color-blind friendly #0072B2 (0,114,178) blue #A2142F (162,20,47) red                                                                                                               
    cat ${HMAWK} | ${hm_script} -vStrokeC="lightgrey" -vStrokeSC="black" -vL="Enrichment(-log10)" -vD="space" -vLegendNum="6" -vC="0:(255,255,255);$MAX:(205,0,0)" -vFontSize=8 > ${OUTFN}.svg
    # convert -density 300 ${OUTFN}.svg ${OUTFN}.png
done

###
## Step 4 (optional): combine heatmaps
###
process_entries() {
    local spacer_counter=1
    local spacer_prefix="Vspacer"
    local spacer_value=5
    if [[ "$1" =~ ^[0-9]+$ ]]; then
        spacer_value=$1
        shift
    fi
    for entry in "$@"; do
        local file=""
        if [[ "$entry" == *:* ]]; then
            file="${entry%%:*}"
            local suffix="${entry#*:}"
            awk -v sfx="$suffix" '{ sub(/\|\|/, "_" sfx "||", $2); print }' "$file"
        else
            file="$entry"
            cat "$file"
        fi
        # Avoid duplicate spacer IDs                                                                                                                                                        
        while grep -q "|-.*${spacer_prefix}${spacer_counter}" "$file"; do
            spacer_counter=$((spacer_counter + 1))
        done
        # Print spacer line after each file except the last                                                                                                                                 
        if [ "$entry" != "${@: -1}" ]; then
            echo -e "|-\t|Spacer|${spacer_value}|${spacer_prefix}${spacer_counter}"
            spacer_counter=$((spacer_counter + 1))
        fi
    done
}

lda_result_dir=${result_dir}
GO_dir=${lda_result_dir}/ldaGOenrichment/GOenrichment
cd $GO_dir # chances are there are : in the file path
process_entries 5 $(for i in $(seq 1 $((${lda_k}))); do echo -n "go_enrichAnalysis_topic${i}_details.hmawk:topic${i} "; done) \
> combined.hmawk
cd -

lda_result_dir=${result_dir}
GO_dir=${lda_result_dir}/ldaGOenrichment/GOenrichment
HMAWK=${GO_dir}/combined.hmawk
OUTFN="${HMAWK%.*}"
MAX=10
cat ${HMAWK} | ${hm_script} -vStrokeC="lightgrey" -vStrokeSC="black" -vL="Enrichment(-log10)" -vD="space" -vLegendNum="6" -vC="0:(255,255,255);$MAX:(205,0,0)" -vFontSize=8 > ${OUTFN}.svg
# convert -density 300 ${OUTFN}.svg ${OUTFN}.png # too large
