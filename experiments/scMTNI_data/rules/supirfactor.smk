

rule run_supirfactor_shallow_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
    benchmark:
        "benchmarks/supirfactor_shallow/{dataset}/{c}.benchmark.txt"
    threads: 10 # to be comparable to other algos
    conda: 
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/supirfactor/env"
    shell:
        """
        python {config[scriptdir]}/SupirFactor/SupirFactor.py --in_exp {input.d} \
        --tf_lst {input.r} \
        --prior_path {input.p} \
        --out_A {output.out_A} \
        --out_P {output.out_P} \
        --model_used Shallow
        """

rule FBS_supirfactor_shallow_get_summary_stat:
    input:
        A=expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --prior_A {input.prior} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_thresholding_inferred_P # signed TFA
        """

rule BH_supirfactor_shallow_get_summary_stat:
    input:
        A=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --prior_A {input.prior} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_thresholding_inferred_P # signed TFA
        """

rule FBS_supirfactor_shallow_get_AUPR:
    input:
        A=expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_shallow/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_thresholding_inferred_P # signed TFA
        """

rule BH_supirfactor_shallow_get_AUPR:
    input:
        A=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_shallow/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_thresholding_inferred_P # signed TFA
        """

# rule supirfactor_shallow_get_clustering_accurarcy:
#     input: 
#         P=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
#         gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
#     output:
#         outstat=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/Seurat"
#     shell:
#         """
#         Rscript {config[scriptdir]}/downstream_analysis/Seurat_clustering.R \
#         {input.P} \
#         {input.gt_celltype} \
#         {config[tflist]} \
#         $'\\t' \
#         F \
#         {output.outstat}
#         """

####
## Hierarchical
####

rule run_supirfactor_hierarchical_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
    benchmark:
        "benchmarks/supirfactor_hierarchical/{dataset}/{c}.benchmark.txt"
    threads: 10 # to be comparable to other algos
    conda: 
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/supirfactor/env"
    shell:
        """
        python {config[scriptdir]}/SupirFactor/SupirFactor.py --in_exp {input.d} \
        --tf_lst {input.r} \
        --prior_path {input.p} \
        --out_A {output.out_A} \
        --out_P {output.out_P} \
        --model_used Hierarchical
        """

rule FBS_supirfactor_hierarchical_get_summary_stat:
    input:
        A=expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --prior_A {input.prior} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P 
        """

rule BH_supirfactor_hierarchical_get_summary_stat:
    input:
        A=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --prior_A {input.prior} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P 
        """

rule FBS_supirfactor_hierarchical_get_AUPR:
    input:
        A=expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_hierarchical/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P 
        """

rule BH_supirfactor_hierarchical_get_AUPR:
    input:
        A=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_hierarchical/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_scaling_inferred_A \
        --not_scaling_inferred_P 
        """

# rule supirfactor_hierarchical_get_clustering_accurarcy:
#     input: 
#         P=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
#         gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
#     output:
#         outstat=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/Seurat"
#     shell:
#         """
#         Rscript {config[scriptdir]}/downstream_analysis/Seurat_clustering.R \
#         {input.P} \
#         {input.gt_celltype} \
#         {config[tflist]} \
#         $'\\t' \
#         F \
#         {output.outstat}
#         """