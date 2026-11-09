
# TODO: change the method of launching matlab. Below we did "module load matlab-r2023a"

rule run_NINCA_lp_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/NINCA_lp/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/NINCA_lp/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/NINCA_lp/{dataset}/{c}.benchmark.txt"
    shell:
        """
        module load matlab-r2022b
        cd {config[scriptdir]}/NINCA
        matlab -nodisplay -nodesktop -r "X_path='{input.d}'; prior_path='{input.p}'; tf_path='{input.r}'; A_out='{output.out_A}'; \
        P_out='{output.out_P}'; run_ncanec_lp_sub_nocvx_colwise(X_path, prior_path, tf_path, A_out, P_out); quit"
        """

rule run_NINCA_quad_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
    benchmark:
        "benchmarks/NINCA_quad/{dataset}/{c}.benchmark.txt"
    shell:
        """
        module load matlab-r2022b
        cd {config[scriptdir]}/NINCA
        matlab -nodisplay -nodesktop -r "X_path='{input.d}'; prior_path='{input.p}'; tf_path='{input.r}'; A_out='{output.out_A}'; \
        P_out='{output.out_P}'; run_ncanec_quad_sub_nocvx_colwise(X_path, prior_path, tf_path, A_out, P_out); quit"
        """

rule FBS_NINCA_quad_get_summary_stat:
    input:
        A=expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell: # TODO: add true_P
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
        --transpose_inferred_A
        """

rule FBS_NINCA_quad_get_AUPR:
    input:
        A=expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NINCA_quad/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
    shell: # TODO: add true_P
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A
        """


rule BH_NINCA_quad_get_summary_stat:
    input:
        A=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --transpose_inferred_A
        """

rule BH_NINCA_quad_get_AUPR:
    input:
        A=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NINCA_quad/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --transpose_inferred_A
        """

# rule NINCA_quad_get_clustering_accurarcy:
#     input: 
#         P=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
#         gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
#     output:
#         outstat=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/Seurat"
#     shell:
#         """
#         Rscript {config[scriptdir]}/downstream_analysis/Seurat_clustering.R \
#         {input.P} \
#         {input.gt_celltype} \
#         {config[tflist]} \
#         ' ' \
#         F \
#         {output.outstat}
#         """