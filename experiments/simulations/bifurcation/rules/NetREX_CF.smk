

rule run_NetREX_CF_TFA_inference:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        out_c=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/confidence.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/NetREX_CF/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}.benchmark.txt"
    threads: 10 # to be comparable to other algos
    conda: # BITFAM and NetREX_CF can have the same conda environment
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/NetREX_CF"
    shell:
        """
        python {config[scriptdir]}/NetREX_CF/NetREX_CF.py --in_exp {input.d} \
        --tf_lst {input.r} \
        --prior_path {input.p} \
        --out_A {output.out_A} \
        --out_P {output.out_P} \
        --out_confidence {output.out_c}
        """

rule NetREX_CF_get_summary_stat:
    input:
        A=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """

rule NetREX_CF_get_clustering_accurarcy:
    input: 
        P=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/P.txt", outdir=config["outdir"]),
        gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
    output:
        outstat=expand("{outdir}/NetREX_CF/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
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