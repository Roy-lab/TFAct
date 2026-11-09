

rule run_NetREX_CF_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        out_c=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}/confidence.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/NetREX_CF/{dataset}/{c}.benchmark.txt"
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

rule run_NetREX_CF_TFA_other_param_inference:
    input:
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        out_c=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}/confidence.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/NetREX_CF_{weight}/{dataset}/{c}.benchmark.txt"
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
        --out_confidence {output.out_c} \
        --weights {wildcards.weight}
        """

rule FBS_NetREX_CF_get_summary_stat:
    input:
        A=expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """


rule BH_NetREX_CF_get_summary_stat:
    input:
        A=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """

rule FBS_NetREX_CF_get_AUPR:
    input:
        A=expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """


rule BH_NetREX_CF_get_AUPR:
    input:
        A=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """


rule FBS_NetREX_CF_other_param_get_summary_stat:
    input:
        A=expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """


rule BH_NetREX_CF_other_param_get_summary_stat:
    input:
        A=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/TIGER/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """

rule FBS_NetREX_CF_other_param_get_AUPR:
    input:
        A=expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF_{{weight}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """


rule BH_NetREX_CF_other_param_get_AUPR:
    input:
        A=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/NetREX_CF_{{weight}}/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"])
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
        --not_thresholding_inferred_P
        """