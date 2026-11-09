
# TODO: change the method of launching matlab. Below we did "module load matlab-r2023a"

rule run_FastNCA_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/FastNCA/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/FastNCA/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/FastNCA/{dataset}/{c}.benchmark.txt"
    threads: 10
    shell:
        """
        module load matlab-r2023a
        cd {config[scriptdir]}/FastNCA
        matlab -nodisplay -nodesktop -r "X_path='{input.d}'; prior_path='{input.p}'; tf_path='{input.r}'; A_out='{output.out_A}'; \
        P_out='{output.out_P}'; run_FastNCA(X_path, prior_path, tf_path, A_out, P_out); quit"
        """

## TODO: we don't have ground truth P
# rule FastNCA_get_summary_stat:
#     input:
#         A=expand("{outdir}/FastNCA/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
#         P=expand("{outdir}/FastNCA/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
#         prior=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
#         true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
#         g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
#     output:
#         expand("{outdir}/FastNCA/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
#         expand("{outdir}/FastNCA/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
#     params:
#         outdir=expand("{outdir}/FastNCA/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}", outdir=config["outdir"])
#     shell:
#         """
#         python {config[scriptdir]}/compare_networks_and_TFA.py \
#         --outdir {params.outdir} \
#         --inferred_A {input.A}  \
#         --inferred_P {input.P} \
#         --true_A {config[genedir]}/{wildcards.edgeprob}/{wildcards.pre}/nw.tsv \
#         --prior_A {input.prior} \
#         --true_P {input.true_P} \
#         --tf_lst {config[tflist]} \
#         --gene {input.g} \
#         --top {wildcards.top_e} \
#         --transpose_inferred_A
#         """

## TODO: we don't have ground truth clustering assignment
# rule FascNCA_get_clustering_accurarcy:
#     input: 
#         P=expand("{outdir}/FastNCA/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
#         gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
#     output:
#         outstat=expand("{outdir}/FastNCA/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
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