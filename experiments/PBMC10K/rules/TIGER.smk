

rule run_TIGER_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/all_regulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/normalized_{{c}}_exp.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/pbmc10k_prior.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/TIGER/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/TIGER/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/TIGER/{dataset}/{c}.benchmark.txt"
    threads: 10 # to be comparable to other algos
    conda: # BITFAM and TIGER can have the same conda environment
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/BITFAM"
    shell:
        """
        Rscript {config[scriptdir]}/TIGER/TIGER.R {input.d} {input.r} {input.p} {output.out_A} {output.out_P}
        """

# rule FBS_TIGER_get_summary_stat:
#     input:
#         A=expand("{outdir}/TIGER/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
#         P=expand("{outdir}/TIGER/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
#         prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
#         gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
#         r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
#         # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
#     output:
#         # expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
#         expand("{outdir}/TIGER/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
#     params:
#         outdir=expand("{outdir}/TIGER/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
#     shell:
#         """
#         python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
#         --outdir {params.outdir} \
#         --inferred_A {input.A}  \
#         --inferred_P {input.P} \
#         --true_A {input.gs} \
#         --prior_A {input.prior} \
#         --tf_lst {input.r} \
#         --gene {input.g} \
#         --top {wildcards.top_e} \
#         --transpose_inferred_A \
#         --not_scaling_inferred_A \
#         --not_scaling_inferred_P
#         """

# rule BH_TIGER_get_summary_stat:
#     input:
#         A=expand("{outdir}/TIGER/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
#         P=expand("{outdir}/TIGER/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
#         prior=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
#         gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
#         r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
#         # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
#     output:
#         # expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
#         expand("{outdir}/TIGER/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
#     params:
#         outdir=expand("{outdir}/TIGER/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
#     shell:
#         """
#         python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
#         --outdir {params.outdir} \
#         --inferred_A {input.A}  \
#         --inferred_P {input.P} \
#         --true_A {input.gs} \
#         --prior_A {input.prior} \
#         --tf_lst {input.r} \
#         --gene {input.g} \
#         --top {wildcards.top_e} \
#         --transpose_inferred_A \
#         --not_scaling_inferred_A \
#         --not_scaling_inferred_P
#         """

# rule TIGER_get_clustering_accurarcy:
#     input: 
#         P=expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
#         gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
#     output:
#         outstat=expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
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