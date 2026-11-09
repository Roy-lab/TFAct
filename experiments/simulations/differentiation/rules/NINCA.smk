
# TODO: change the method of launching matlab. Below we did "module load matlab-r2023a"

rule run_NINCA_lp_TFA_inference:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/NINCA_lp/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/NINCA_lp/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/NINCA_lp/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}.benchmark.txt"
    shell:
        """
        module load matlab-r2022b
        cd {config[scriptdir]}/NINCA
        matlab -nodisplay -nodesktop -r "X_path='{input.d}'; prior_path='{input.p}'; tf_path='{input.r}'; A_out='{output.out_A}'; \
        P_out='{output.out_P}'; run_ncanec_lp_sub_nocvx_colwise(X_path, prior_path, tf_path, A_out, P_out); quit"
        """

rule run_NINCA_quad_TFA_inference:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/NINCA_quad/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}.benchmark.txt"
    shell:
        """
        module load matlab-r2022b
        cd {config[scriptdir]}/NINCA
        matlab -nodisplay -nodesktop -r "X_path='{input.d}'; prior_path='{input.p}'; tf_path='{input.r}'; A_out='{output.out_A}'; \
        P_out='{output.out_P}'; run_ncanec_quad_sub_nocvx_colwise(X_path, prior_path, tf_path, A_out, P_out); quit"
        """

rule NINCA_quad_get_summary_stat:
    input:
        A=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/compare_networks_and_TFA.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {config[genedir]}/{wildcards.edgeprob}/{wildcards.pre}/nw.tsv \
        --prior_A {input.prior} \
        --true_P {input.true_P} \
        --tf_lst {config[tflist]} \
        --gene {input.g} \
        --top {wildcards.top_e} \
        --transpose_inferred_A \
        --not_thresholding_inferred_P # can have negative signs
        """

rule NINCA_quad_get_clustering_accurarcy:
    input: 
        P=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
    output:
        outstat=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/Seurat"
    shell:
        """
        Rscript {config[scriptdir]}/downstream_analysis/Seurat_clustering.R \
        {input.P} \
        {input.gt_celltype} \
        {config[tflist]} \
        ' ' \
        F \
        {output.outstat}
        """

rule NINCA_quad_get_AUPR:
    input:
        A=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gs=expand("{genedir}/{{edgeprob}}/{{pre}}/nw.tsv", genedir=config["genedir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR_TF.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/", outdir=config["outdir"]),
    shell: 
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --transpose_inferred_A \
        --not_outputing_top_edges
        """

rule NINCA_quad_get_AUPR_prior:
    input:
        A=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR_TF.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NINCA_quad/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/", outdir=config["outdir"]),
    shell: 
        """
        python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
        --outdir {params.outdir} \
        --inferred_A {input.A}  \
        --inferred_P {input.P} \
        --true_A {input.gs} \
        --tf_lst {input.r} \
        --gene {input.g} \
        --transpose_inferred_A \
        --not_outputing_top_edges
        """