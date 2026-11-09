

rule run_supirfactor_shallow_TFA_inference:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/supirfactor_shallow/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}.benchmark.txt"
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

rule supirfactor_shallow_get_summary_stat:
    input:
        A=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_thresholding_inferred_P # signed TFA
        """

rule supirfactor_shallow_get_clustering_accurarcy:
    input: 
        P=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
    output:
        outstat=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/Seurat"
    shell:
        """
        Rscript {config[scriptdir]}/downstream_analysis/Seurat_clustering.R \
        {input.P} \
        {input.gt_celltype} \
        {config[tflist]} \
        $'\\t' \
        F \
        {output.outstat}
        """

rule supirfactor_shallow_get_AUPR:
    input:
        A=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gs=expand("{genedir}/{{edgeprob}}/{{pre}}/nw.tsv", genedir=config["genedir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR_TF.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/", outdir=config["outdir"]),
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
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_thresholding_inferred_P \
        --not_outputing_top_edges
        """

rule supirfactor_shallow_get_AUPR_prior:
    input:
        A=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR_TF.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_shallow/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/", outdir=config["outdir"]),
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
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_thresholding_inferred_P \
        --not_outputing_top_edges
        """

####
## Hierarchical
####

rule run_supirfactor_hierarchical_TFA_inference:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/supirfactor_hierarchical/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}.benchmark.txt"
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

rule supirfactor_hierarchical_get_summary_stat:
    input:
        A=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_scaling_inferred_A \
        --not_scaling_inferred_P 
        """


rule supirfactor_hierarchical_get_clustering_accurarcy:
    input: 
        P=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
    output:
        outstat=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/Seurat"
    shell:
        """
        Rscript {config[scriptdir]}/downstream_analysis/Seurat_clustering.R \
        {input.P} \
        {input.gt_celltype} \
        {config[tflist]} \
        $'\\t' \
        F \
        {output.outstat}
        """

rule supirfactor_hierarchical_get_AUPR:
    input:
        A=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gs=expand("{genedir}/{{edgeprob}}/{{pre}}/nw.tsv", genedir=config["genedir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR_TF.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/", outdir=config["outdir"]),
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
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_outputing_top_edges
        """

rule supirfactor_hierarchical_get_AUPR_prior:
    input:
        A=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR_TF.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/supirfactor_hierarchical/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/", outdir=config["outdir"]),
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
        --not_scaling_inferred_A \
        --not_scaling_inferred_P \
        --not_outputing_top_edges
        """