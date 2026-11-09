###
# Baseline: merlin without TFA
###
rule merlin_create_expression_only_datasets:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
    output:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_replaced/no_TFA/{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
    params:
        indir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}", datadir=config["datadir"]),
        outdir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_replaced/no_TFA", datadir=config["datadir"]),
        num_partitions=run_num_subsamps,
    shell:
        """
        for i in $(seq 0 $(({params.num_partitions}-1)) ); do
            python {config[scriptdir]}/downstream_network_inference/create_merlin_datasets_with_TFA_features.py \
                --exp {params.indir}/dataset.${{i}}.txt \
                --expression_only \
                --out {params.outdir}/${{i}}.txt
        done
        """


###
# In the first scenario, we created new datasets by replacing expression with TFA
###

rule merlin_create_datasets_w_exp_replaced_by_tfa:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
        tfa=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA/{{algo}}/{subsamp}.tsv.gz", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_replaced/{{algo}}/{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
    params:
        indir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}", datadir=config["datadir"]),
        outdir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_replaced/{{algo}}", datadir=config["datadir"]),
        num_partitions=run_num_subsamps,
    shell:
        """
        for i in $(seq 0 $(({params.num_partitions}-1)) ); do
            python {config[scriptdir]}/downstream_network_inference/create_merlin_datasets_with_TFA_features.py \
                --exp {params.indir}/dataset.${{i}}.txt \
                --tf_lst {input.r} \
                --out {params.outdir}/${{i}}.txt \
                --tfa {params.indir}/TFA/{wildcards.algo}/${{i}}.tsv.gz \
                --replace 
        done
        """

rule merlin_create_regularized_nca_datasets_w_exp_replaced_by_tfa:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
        tfa=expand("{TFAdir}/regularized_nca/{{dataset}}/{{c}}_stability/P.txt", \
            TFAdir=config["TFAdir"], subsamp=str(run_num_subsamps-1)),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_replaced/regularized_nca/{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
    params:
        indir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}", datadir=config["datadir"]),
        outdir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_replaced/regularized_nca", datadir=config["datadir"]),
        num_partitions=run_num_subsamps,
    shell:
        """
        for i in $(seq 0 $(({params.num_partitions}-1)) ); do
            python {config[scriptdir]}/downstream_network_inference/create_merlin_datasets_with_TFA_features.py \
                --exp {params.indir}/dataset.${{i}}.txt \
                --tf_lst {input.r} \
                --out {params.outdir}/${{i}}.txt \
                --tfa {input.tfa} \
                --replace \
                --tfa_transposed_no_header
        done
        """

###
# In the second scenario, we created new datasets by concatenating TFA values
###


rule merlin_create_datasets_w_tfa_concatenated:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
        tfa=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA/{{algo}}/{subsamp}.tsv.gz", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_concat/{{algo}}/{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
    params:
        indir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}", datadir=config["datadir"]),
        outdir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_concat/{{algo}}", datadir=config["datadir"]),
        num_partitions=run_num_subsamps,
    shell:
        """
        for i in $(seq 0 $(({params.num_partitions}-1)) ); do
            python {config[scriptdir]}/downstream_network_inference/create_merlin_datasets_with_TFA_features.py \
                --exp {params.indir}/dataset.${{i}}.txt \
                --tf_lst {input.r} \
                --out {params.outdir}/${{i}}.txt \
                --tfa {params.indir}/TFA/{wildcards.algo}/${{i}}.tsv.gz
        done
        """

rule merlin_create_regularized_nca_datasets_w_tfa_concat:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
        tfa=expand("{TFAdir}/regularized_nca/{{dataset}}/{{c}}_stability/P.txt", \
            TFAdir=config["TFAdir"], subsamp=str(run_num_subsamps-1)),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_concat/regularized_nca/{subsamp}.txt", \
            datadir=config["datadir"], subsamp=str(run_num_subsamps-1)),
    params:
        indir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}", datadir=config["datadir"]),
        outdir=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_concat/regularized_nca", datadir=config["datadir"]),
        num_partitions=run_num_subsamps,
    shell:
        """
        for i in $(seq 0 $(({params.num_partitions}-1)) ); do
            python {config[scriptdir]}/downstream_network_inference/create_merlin_datasets_with_TFA_features.py \
                --exp {params.indir}/dataset.${{i}}.txt \
                --tf_lst {input.r} \
                --out {params.outdir}/${{i}}.txt \
                --tfa {input.tfa} \
                --tfa_transposed_no_header
        done
        """

rule merlin_create_TFA_prior:
    input:
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network_TFA.txt", datadir=config["datadir"]),
    shell:
        """
        awk '{{OFS="\\t"}}{{print $1 "_TFA", $2, $3}}1' {input.p} > {output.p}
        """

rule merlin_create_prior_config:
    input:
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network_TFA.txt", datadir=config["datadir"]),
    output:
        p_config=expand("{datadir}/{{dataset}}/inputdata/configdir/{{c}}_{{prior_conf}}.txt", datadir=config["datadir"]),
    shell:
        """
        echo -e "{wildcards.dataset}_{wildcards.c}\\t{input.p}\\t{wildcards.prior_conf}" > {output.p_config}
        """

# actually this gives a list of gene expressed in at least merlin_regulator_min_cell cells in cluster c.
rule merlin_filter_regulator:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
    output:
        out=expand("{datadir}/{{dataset}}/inputdata/merlin/{{c}}_expressed_genes.txt", datadir=config["datadir"]),
    params:
        min_cell=merlin_regulator_min_cell
    shell:
        """
        python {config[scriptdir]}/downstream_network_inference/merlin_get_expressed_genes.py --exp {input.d} \
            --out {output.out} --thr {params.min_cell}
        """

rule merlin_FBS_combine_regulator:
    input:
        expand("{datadir}/processed_FBS_reprogramming/inputdata/merlin/{c}_expressed_genes.txt", datadir=config["datadir"], c = FBS_clusters),
    output:
        expand("{datadir}/processed_FBS_reprogramming/inputdata/merlin/expressed_genes_combined.txt", datadir=config["datadir"]),
    shell:
        """
        cat {input} | sort -u > {output}
        """

rule merlin_BH_combine_regulator:
    input:
        expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/merlin/{c}_expressed_genes.txt", datadir=config["datadir"], c = BH_clusters),
    output:
        expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/merlin/expressed_genes_combined.txt", datadir=config["datadir"]),
    shell:
        """
        cat {input} | sort -u > {output}
        """

rule merlin_RFH_combine_regulator:
    input:
        expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/merlin/{c}_expressed_genes.txt", datadir=config["datadir"], c = RFH_fine_clusters),
    output:
        expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/merlin/expressed_genes_combined.txt", datadir=config["datadir"]),
    shell:
        """
        cat {input} | sort -u > {output}
        """

# Here, we include all regulators, not just those have TFA values or in the prior networks.
# if a regulator does not have TFA values, only include regulators that appear in at least 50 cells in some cell cluster
rule merlin_create_TFA_regulators_FBS:
    input:
        exp_r=expand("{datadir}/FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        exp_g=expand("{datadir}/processed_FBS_reprogramming/inputdata/merlin/expressed_genes_combined.txt", datadir=config["datadir"]),
        tfa_r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
    shell:
        """
        comm -12 <(sort {input.exp_g}) <(awk -v cluster="_{wildcards.c}" '{{sub(cluster "$", "", $1); print $1}}' {input.exp_r} | sort ) > {input.exp_r}.tmp
        awk 'FILENAME == ARGV[1] {{print;next}} FILENAME == ARGV[2] {{print $1 "_TFA"}}1' {input.exp_r}.tmp {input.tfa_r} | sort -u > {output}
        rm {input.exp_r}.tmp
        """

rule merlin_create_TFA_regulators_BH:
    input:
        exp_r=expand("{datadir}/Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        exp_g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/merlin/expressed_genes_combined.txt", datadir=config["datadir"]),
        tfa_r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
    shell:
        """
        comm -12 <(sort {input.exp_g}) <(awk -v cluster="_{wildcards.c}" '{{sub(cluster "$", "", $1); print $1}}' {input.exp_r} | sort ) > {input.exp_r}.tmp
        awk 'FILENAME == ARGV[1] {{print;next}} FILENAME == ARGV[2] {{print $1 "_TFA"}}1' {input.exp_r}.tmp {input.tfa_r} | sort -u > {output}
        rm {input.exp_r}.tmp
        """

rule merlin_create_TFA_regulators_RFH:
    input:
        exp_r=expand("{datadir}/Ranzoni_fetal_hematopoiesis/fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        exp_g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/merlin/expressed_genes_combined.txt", datadir=config["datadir"]),
        tfa_r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
    shell:
        """
        comm -12 <(sort {input.exp_g}) <(awk -v cluster="_{wildcards.c}" '{{sub(cluster "$", "", $1); print $1}}' {input.exp_r} | sort ) > {input.exp_r}.tmp
        awk 'FILENAME == ARGV[1] {{print;next}} FILENAME == ARGV[2] {{print $1 "_TFA"}}1' {input.exp_r}.tmp {input.tfa_r} | sort -u > {output}
        rm {input.exp_r}.tmp
        """

##
## Note, as the current snakemake condor profile does not enable flocking, we run merlin directly on condor instead.
## the rules assume that condor jobs have been run separately, and the results have been organized using the organize_merlin_condor_output.sh script.
##
rule merlin_make_consensus_file:
    input:
        outnet=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/{subsamp}/fold0/prediction_k300.txt", outdir=config["outdir"], subsamp=str(run_num_subsamps-1)),
    output:
        outnet_set=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/stability/network_files.txt", outdir=config["outdir"]),
        out_consensus=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]),
        consensus_dist=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/stability/consensus_dist.txt", outdir=config["outdir"]),
    params:
        outnet_format=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/*/fold0/prediction_k300.txt", outdir=config["outdir"]),
        estEdge=estEdge,
        out_consensus_prefix=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/stability/consensus_", outdir=config["outdir"])
    shell:
        """
        ls {params.outnet_format} > {output.outnet_set} 
        {params.estEdge} {output.outnet_set} 0 {params.out_consensus_prefix} alledges
        """


##
# Evaluation
# For TFA concat: if we have an edge TF->target and an edge TF_TFA->target, then we construct the inferred networks in 4 ways: 1. remove TFA edges. 2. remove edges without TFA. 3. take the max confidence. 4. take the min confidence
# We also have TFA->TFA edges, and TF->TFA edges
##
rule merlin_FBS_replaced_get_AUPR:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]), 
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --A_edge_lst
        """

## 
rule merlin_FBS_concat_get_AUPR:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]), 
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --A_edge_lst \
        --separate_eval_by_edge_types
        """


rule merlin_BH_replaced_get_AUPR:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]), 
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --A_edge_lst
        """

## 
rule merlin_BH_concat_get_AUPR:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]), 
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --A_edge_lst \
        --separate_eval_by_edge_types
        """

rule merlin_RFH_replaced_get_AUPR:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]), 
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --A_edge_lst
        """

## 
rule merlin_RFH_concat_get_AUPR:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]), 
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --A_edge_lst \
        --separate_eval_by_edge_types
        """

# rule merlin_BH_get_AUPR:
#     input:
#         A=expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/network.txt", outdir=config["outdir"]),
#         gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
#         r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
#     output:
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
#     params:
#         outdir=expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
#     shell: # TODO add true_P;
#         """
#         python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
#         --outdir {params.outdir} \
#         --inferred_A {input.A}  \
#         --true_A {input.gs} \
#         --tf_lst {input.r} \
#         --gene {input.g} \
#         --top {wildcards.top_e} \
#         --A_edge_lst
#         """

# rule merlin_RFH_get_AUPR:
#     input:
#         A=expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/network.txt", outdir=config["outdir"]),
#         gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
#         r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
#     output:
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
#     params:
#         outdir=expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
#     shell: # TODO add true_P;
#         """
#         python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
#         --outdir {params.outdir} \
#         --inferred_A {input.A}  \
#         --true_A {input.gs} \
#         --tf_lst {input.r} \
#         --gene {input.g} \
#         --top {wildcards.top_e} \
#         --A_edge_lst
#         """

rule merlin_FBS_replaced_get_F:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(merlin_algos + ["no_TFA"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        """

rule merlin_FBS_concat_get_F:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(merlin_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst \
            --separate_eval_by_edge_types
        """

rule merlin_BH_replaced_get_F:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(merlin_algos + ["no_TFA"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        """

rule merlin_BH_concat_get_F:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(merlin_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst \
            --separate_eval_by_edge_types
        """

rule merlin_RFH_replaced_get_F:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(merlin_algos + ["no_TFA"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_replaced/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        """

rule merlin_RFH_concat_get_F:
    input:
        A=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/consensus_alledge.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/max_conf_stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(merlin_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/merlin_{{algo}}_concat/c{{prior_conf}}_p{{sparsity}}_r{{module}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst \
            --separate_eval_by_edge_types
        """

# rule merlin_BH_get_F:
#     input:
#         A=expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/network.txt", outdir=config["outdir"]),
#         prior=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
#         gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
#         r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
#     output:
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
#     wildcard_constraints:
#         algo="|".join(merlin_algos)
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
#     params:
#         outdir=expand("{outdir}/algos/merlin_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
#     shell:
#         """
#         python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
#             --outdir {params.outdir} \
#             --inferred_A {input.A}  \
#             --true_A {input.gs} \
#             --prior_A {input.prior} \
#             --tf_lst {input.r} \
#             --gene {input.g} \
#             --top {wildcards.top_e} \
#             --A_edge_lst
#         """

# rule merlin_RFH_get_F:
#     input:
#         A=expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/network.txt", outdir=config["outdir"]),
#         prior=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
#         gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
#         r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
#     output:
#         expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
#     wildcard_constraints:
#         algo="|".join(merlin_algos)
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
#     params:
#         outdir=expand("{outdir}/algos/merlin_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
#     shell:
#         """
#         python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
#             --outdir {params.outdir} \
#             --inferred_A {input.A}  \
#             --true_A {input.gs} \
#             --prior_A {input.prior} \
#             --tf_lst {input.r} \
#             --gene {input.g} \
#             --top {wildcards.top_e} \
#             --A_edge_lst
#         """

## the following is deprecated
# rule run_merlin_p:
#     input:
#         input_ds=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA_algo_{{ds_action}}/{{algo}}/{{subsamp}}.txt", \
#             datadir=config["datadir"]),
#         prior_config=expand("{datadir}/{{dataset}}/inputdata/configdir/{{c}}_{{prior_conf}}.txt", datadir=config["datadir"]),
#         reg=expand("{datadir}/{{dataset}}/inputdata/{{c}}_complete_allregulators_TFA.txt", datadir=config["datadir"]),
#     output:
#         outnet=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c_{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/{{subsamp}}/fold0/prediction_k300.txt", outdir=config["outdir"]),
#         outmod=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c_{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/{{subsamp}}/fold0/modules.txt", outdir=config["outdir"]),
#     threads: 1
#     resources:
#         mem_mb=20000,
#         disk_mb=5000
#     params:
#         code_path="/mnt/dv/wid/projects5/Roy-singlecell/netinf_algo_comparison/bin/MERLIN/merlin-p_self_checkpoint/merlin", # change this
#         LD_LIB="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/envs/IBPG/lib", # change this
#         outdir=expand("{outdir}/algos/merlin_{{algo}}_{{ds_action}}/c_{{prior_conf}}_p{{sparsity}}_r{{module}}/{{dataset}}/{{c}}/{{subsamp}}", outdir=config["outdir"]),
#     shell:
#         """
#         export LD_LIBRARY_PATH=${{LD_LIBRARY_PATH}}:{params.LD_LIB}
#         {params.code_path} \
#             -d {input.input_ds} \
#             -l {input.reg} \
#             -o {params.outdir} \
#             -q {input.prior_config}	\
#             -p {wildcards.sparsity} \
#             -r {wildcards.module}
#         """