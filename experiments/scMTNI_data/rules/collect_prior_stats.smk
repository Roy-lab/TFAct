rule FBS_prior_get_AUPR:
    input:
        A=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/prior/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/prior/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/prior/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/prior/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell: # TODO add true_P;
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

rule BH_prior_get_AUPR:
    input:
        A=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/prior/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/prior/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/prior/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/prior/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell: # TODO add true_P;
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

rule RFH_prior_get_AUPR:
    input:
        A=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/prior/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/prior/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/prior/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/prior/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
    shell: # TODO add true_P;
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