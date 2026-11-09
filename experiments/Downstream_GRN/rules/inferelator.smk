rule run_inferelator_plain:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{{subsamp}}.txt", \
            datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/inferelator.{{c}}_network.tsv.gz", datadir=config["datadir"]),
    output:
        grn=expand("{outdir}/algos/inferelator_default/{{dataset}}/{{c}}/{{subsamp}}/network.tsv.gz", outdir=config["outdir"]),
    params:
        outdir=expand("{outdir}/algos/inferelator_default/{{dataset}}/{{c}}/{{subsamp}}", outdir=config["outdir"]),
    benchmark:
        "benchmarks/inferelator_default/{dataset}.{c}.{subsamp}.benchmark.txt"
    conda: "/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/env/inferelator3"
    threads: 5
    resources:
        mem_mb=5000,
        disk_mb=10000
    shell:
        """
        python {config[scriptdir]}/downstream_network_inference/run_inferelator3.py \
            --exp {input.d} \
            --tf_lst {input.r} \
            --prior_A {input.p} \
            --gs_A {input.p} \
            --outdir {params.outdir} \
            --threads {threads}
        tar czf {params.outdir}/inferelator_model.h5ad.tar.gz {params.outdir}/inferelator_model.h5ad
        rm {params.outdir}/inferelator_model.h5ad
        """

rule run_inferelator_no_TFA:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{{subsamp}}.txt", \
            datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/inferelator.{{c}}_network.tsv.gz", datadir=config["datadir"]),
    output:
        grn=expand("{outdir}/algos/inferelator_no_TFA/{{dataset}}/{{c}}/{{subsamp}}/network.tsv.gz", outdir=config["outdir"]),
    params:
        outdir=expand("{outdir}/algos/inferelator_no_TFA/{{dataset}}/{{c}}/{{subsamp}}", outdir=config["outdir"]),
    benchmark:
        "benchmarks/inferelator_no_TFA/{dataset}.{c}.{subsamp}.benchmark.txt"
    conda: "/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/env/inferelator3"
    threads: 5
    resources:
        mem_mb=5000,
        disk_mb=10000
    shell:
        """
        python {config[scriptdir]}/downstream_network_inference/run_inferelator3.py \
            --exp {input.d} \
            --tf_lst {input.r} \
            --prior_A {input.p} \
            --gs_A {input.p} \
            --outdir {params.outdir} \
            --no_TFA \
            --threads {threads}
        tar czf {params.outdir}/inferelator_model.h5ad.tar.gz {params.outdir}/inferelator_model.h5ad
        rm {params.outdir}/inferelator_model.h5ad
        """

rule run_inferelator_external_tfa:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{{subsamp}}.txt", \
            datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]), # to make algos comparable, we don't use the union of tfs across cell types
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/inferelator.{{c}}_network.tsv.gz", datadir=config["datadir"]),
        tfa=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA/{{algo}}/{{subsamp}}.tsv.gz", \
            datadir=config["datadir"]),
    output:
        grn=expand("{outdir}/algos/inferelator_{{algo}}/{{dataset}}/{{c}}/{{subsamp}}/network.tsv.gz", outdir=config["outdir"]),
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/{{dataset}}/{{c}}/{{subsamp}}", outdir=config["outdir"]),
    benchmark:
        "benchmarks/inferelator_{algo}/{dataset}.{c}.{subsamp}.benchmark.txt"
    wildcard_constraints:
        algo='^(?!regularized_nca_tfa_only$)(?!regularized_nca_tfa_and_exp$).+'
    conda: "/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/env/inferelator3"
    threads: 5
    resources:
        mem_mb=5000,
        disk_mb=10000
    shell:
        """
        python {config[scriptdir]}/downstream_network_inference/run_inferelator3.py \
            --exp {input.d} \
            --tf_lst {input.r} \
            --prior_A {input.p} \
            --gs_A {input.p} \
            --outdir {params.outdir} \
            --tfa {input.tfa} \
            --threads {threads}
        tar czf {params.outdir}/inferelator_model.h5ad.tar.gz {params.outdir}/inferelator_model.h5ad
        rm {params.outdir}/inferelator_model.h5ad
        """

rule run_inferelator_regularized_nca:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/dataset.{{subsamp}}.txt", \
            datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA/{{algo}}/regulators.txt", \
            datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/inferelator.{{c}}_network.tsv.gz", datadir=config["datadir"]),
        tfa=expand("{datadir}/{{dataset}}/inputdata/subsamples/{{c}}/TFA/{{algo}}/{{subsamp}}.tsv.gz", \
            datadir=config["datadir"]),
    output:
        grn=expand("{outdir}/algos/inferelator_{{algo}}/{{dataset}}/{{c}}/{{subsamp}}/network.tsv.gz", outdir=config["outdir"]),
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/{{dataset}}/{{c}}/{{subsamp}}", outdir=config["outdir"]),
    benchmark:
        "benchmarks/inferelator_{algo}/{dataset}.{c}.{subsamp}.benchmark.txt"
    wildcard_constraints:
        algo="|".join(['regularized_nca_tfa_only', 'regularized_nca_tfa_and_exp'])
    conda: "/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA/env/inferelator3"
    threads: 5
    resources:
        mem_mb=5000,
        disk_mb=10000
    shell:
        """
        python {config[scriptdir]}/downstream_network_inference/run_inferelator3.py \
            --exp {input.d} \
            --tf_lst {input.r} \
            --prior_A {input.p} \
            --gs_A {input.p} \
            --outdir {params.outdir} \
            --tfa {input.tfa} \
            --threads {threads}
        tar czf {params.outdir}/inferelator_model.h5ad.tar.gz {params.outdir}/inferelator_model.h5ad
        rm {params.outdir}/inferelator_model.h5ad
        """

rule stability_selection:
    input:
        grn=expand("{outdir}/algos/inferelator_{{algo}}/{{dataset}}/{{c}}/{subsamp}/network.tsv.gz", outdir=config["outdir"], subsamp=[str(_) for _ in range(run_num_subsamps)]),
    output:
        grn_stability=expand("{outdir}/algos/inferelator_{{algo}}/{{dataset}}/{{c}}/stability/network.txt", outdir=config["outdir"]),
    params:
        indir=expand("{outdir}/algos/inferelator_{{algo}}/{{dataset}}/{{c}}", outdir=config["outdir"]),
        num_sub=run_num_subsamps
    shell:
        """
        python {config[scriptdir]}/downstream_network_inference/inferelator_stability_selection.py \
            --indir {params.indir} \
            --num_subsamp {params.num_sub} \
            --outfile {output.grn_stability}
        """

##
# Evaluation
##
rule inferelator_FBS_get_AUPR:
    input:
        A=expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/network.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
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

rule inferelator_BH_get_AUPR:
    input:
        A=expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/network.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
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

rule inferelator_RFH_get_AUPR:
    input:
        A=expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/network.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
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

rule inferelator_FBS_get_F:
    input:
        A=expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/network.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(inferelator_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/processed_FBS_reprogramming/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
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

rule inferelator_BH_get_F:
    input:
        A=expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/network.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(inferelator_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
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

rule inferelator_RFH_get_F:
    input:
        A=expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/network.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(inferelator_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/algos/inferelator_{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}/stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
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