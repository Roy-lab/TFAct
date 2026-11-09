# ref: https://github.com/nyu-dl/pmf-grn
# 
####################
## Step 1: convert expression and prior to h5ad
####################

rule pmf_grn_convert_exp_to_h5ad:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        d=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/{{c,[^/]+}}.h5ad", datadir=config["datadir"])
        # p=expand("{datadir}/{{edgeprob}}/{{ldak}}/pmf_grn/{{pre}}_prior.h5ad", datadir=config["datadir"])
    conda: # TODO
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/pmf-grn/env" # built from the supplied environment yml file. Though the file itself can be put here.
    shell:
        """
        python {config[scriptdir]}/pmf-grn/convert_exp_prior_to_h5ad.py \
        --in_exp {input.d} \
        --in_prior {input.p} \
        --tf_lst {input.r} \
        --out_exp {output.d} 
        """

rule pmf_grn_convert_prior_to_h5ad:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        r=expand("{datadir}/{{dataset}}/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
    output:
        # d=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/pmf_grn/{{exp}}.h5ad", datadir=config["datadir"])
        p=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/priors/{{c,[^/]+}}.h5ad", datadir=config["datadir"])
    conda: # TODO
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/pmf-grn/env" # built from the supplied environment yml file. Though the file itself can be put here.
    shell:
        """
        python {config[scriptdir]}/pmf-grn/convert_exp_prior_to_h5ad.py \
        --in_exp {input.d} \
        --in_prior {input.p} \
        --tf_lst {input.r} \
        --out_prior {output.p}
        """

####################
## Step 2: Split training and validation prior.
####################

rule pmf_grn_split_prior:
    input:
        p=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/priors/{{c}}.h5ad", datadir=config["datadir"])
    output:
        t=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/{{c}}/training_prior_1.h5ad", datadir=config["datadir"]),
        v=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/{{c}}/validation_gs_1.tsv", datadir=config["datadir"])
    params:
        outdir=config["datadir"] + "/{dataset}/inputdata/pmf_grn/{c}"
    conda: # TODO
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/pmf-grn/env" # built from the supplied environment yml file. Though the file itself can be put here.
    shell:
        """
        python {config[scriptdir]}/pmf-grn/split_prior.py \
        --original-prior-path {input.p} \
        --output-dir {params.outdir} \
        --val-frac 0.2 
        """

####################
## Step 3: Hyperparameter search
####################
rule pmf_grn_hyper_search:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/{{c}}.h5ad", datadir=config["datadir"]),
        t=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/{{c}}/training_prior_1.h5ad", datadir=config["datadir"]),
        v=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/{{c}}/validation_gs_1.tsv", datadir=config["datadir"])
    benchmark:
        "benchmarks/pmf_grn_hyper_search/{dataset}/{c}/{fold_num}.txt"
    output:
        best_param=expand("{outdir}/pmf_grn/experiments/hpsearch/{{dataset}}/{{c}}/{{fold_num}}/best_val_auprc_iter.params", outdir=config["outdir"]),
        hyperparam=expand("{outdir}/pmf_grn/experiments/hpsearch/{{dataset}}/{{c}}/{{fold_num}}/hyperparams_used.txt", outdir=config["outdir"])
    threads: 10
    params:
        outdir=config["outdir"] + "/pmf_grn/experiments/hpsearch/{dataset}/{c}/{fold_num}"
    conda: # TODO
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/pmf-grn/env"
    shell: 
        """
        possible_A_stds=({config[algo][pmf_grn][possible_A_stds]})
        A_std=${{possible_A_stds["$((RANDOM % "${{#possible_A_stds[@]}}"))"]}}
        possible_U_stds=({config[algo][pmf_grn][possible_U_stds]})
        U_std=${{possible_U_stds["$((RANDOM % "${{#possible_U_stds[@]}}"))"]}}
        possible_annealing_factors=({config[algo][pmf_grn][possible_annealing_factors]})
        annealing_factor=${{possible_annealing_factors["$((RANDOM % "${{#possible_annealing_factors[@]}}"))"]}}
        echo -e "A_std\tU_std\tannealing_factor" > {output.hyperparam}
        echo -e "${{A_std}}\t${{U_std}}\t${{annealing_factor}}" >> {output.hyperparam}

        cd {config[algo][pmf_grn][dir_path]}

        python -m training.train --expression-path {input.d}  \
        --gene-tf-prior-path {input.t} --gd-type full_dataset \
        --batch-size 0 --num-epochs {config[algo][pmf_grn][num_hyper_search_epoch]} \
        --save-period {config[algo][pmf_grn][num_hyper_search_epoch]} --save-dir {params.outdir} --module-name Gaussian3M --lr 0.1 --min-lr 0.1 --clip-norm 0.0001 \
        --min-prior-hparam 0.0005 --max-prior-hparam 0.9995 --prior-std-logit-A ${{A_std}} --prior-std-log-U ${{U_std}} \
        --val-after 100 --first-100-val-after 10 --val-gold-standard-path {input.v} \
        --num-auprc-logit-sampling-iters 100 \
        --val-kl-anneal-factor 1 --num-annealing-iters 999999 --initial-annealing-factor ${{annealing_factor}} --guide-max-std-log-U 2 
        """

####################
## Step 4: Find the best hyperparameters
####################

rule pmf_grn_find_best_hyper:
    input:
        best_param=expand("{outdir}/pmf_grn/experiments/hpsearch/{{dataset}}/{{c}}/{fold_num}/best_val_auprc_iter.params", \
            fold_num=[_ for _ in range(config["algo"]["pmf_grn"]["num_experiments"])], outdir=config["outdir"])
    output:
        result=expand("{outdir}/pmf_grn/experiments/hpsearch/{{dataset}}/{{c}}/hp_results.txt", outdir=config["outdir"])
    params:
        parent_d=config["outdir"] + "/pmf_grn/experiments/hpsearch/{dataset}/{c}"
    shell:
        """
        python {config[scriptdir]}/pmf-grn/get_hp_search_results.py --parent-dnames {params.parent_d} --output-path {output.result}
        """


####################
## Step 5: Run pmf_grn with full prior and the best hyperparameters
####################

rule run_pmf_grn_TFA_inference:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/{{c}}.h5ad", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/inputdata/pmf_grn/priors/{{c}}.h5ad", datadir=config["datadir"]),
        hp_result=expand("{outdir}/pmf_grn/experiments/hpsearch/{{dataset}}/{{c}}/hp_results.txt", outdir=config["outdir"])
    benchmark:
        "benchmarks/pmf_grn/{dataset}/{c}.benchmark.txt"
    output:
        best_param=expand("{outdir}/pmf_grn/experiments/full_prior/{{dataset}}/{{c}}/best_iter.params", outdir=config["outdir"])
    threads: 10
    params:
        hp_indir=config["outdir"] + "/pmf_grn/experiments/hpsearch/{dataset}/{c}",
        outdir=config["outdir"] + "/pmf_grn/experiments/full_prior/{dataset}/{c}"
    conda: # TODO
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/pmf-grn/env"
    shell: 
        """
        best_fold_num=$(awk -F' ' '$0~/***Best/{{print $NF}}' {input.hp_result})
        best_fold_hp_path={params.hp_indir}/${{best_fold_num}}/hyperparams_used.txt
        A_std=$(awk 'NR==2{{print $1}}' ${{best_fold_hp_path}})
        U_std=$(awk 'NR==2{{print $2}}' ${{best_fold_hp_path}})
        annealing_factor=$(awk 'NR==2{{print $3}}' ${{best_fold_hp_path}})

        cd {config[algo][pmf_grn][dir_path]}

        python -m training.train --expression-path {input.d}  \
        --gene-tf-prior-path {input.p} --gd-type full_dataset \
        --batch-size 0 --num-epochs {config[algo][pmf_grn][num_hyper_search_epoch]} \
        --save-period {config[algo][pmf_grn][num_hyper_search_epoch]} --save-dir {params.outdir} --module-name Gaussian3M --lr 0.1 --min-lr 0.1 --clip-norm 0.0001 \
        --min-prior-hparam 0.0005 --max-prior-hparam 0.9995 --prior-std-logit-A ${{A_std}} --prior-std-log-U ${{U_std}} \
        --val-after 100 --first-100-val-after 10 \
        --num-auprc-logit-sampling-iters 100 --val-loss-type neg_mll --num-mll-samples 10 \
        --val-kl-anneal-factor 1 --num-annealing-iters 999999 --initial-annealing-factor ${{annealing_factor}} --guide-max-std-log-U 2 
        """

####################
## Step 6: Get Consensus networks
####################
rule pmf_grn_get_consensus_net:
    input:
        expand("{outdir}/pmf_grn/experiments/full_prior/{{dataset}}/{{c}}/best_iter.params", outdir=config["outdir"])
    output:
        expand("{outdir}/pmf_grn/{{dataset}}/{{c}}/inferred_consensus_grn.tsv", outdir=config["outdir"]),
        expand("{outdir}/pmf_grn/{{dataset}}/{{c}}/inferred_consensus_tfa.tsv", outdir=config["outdir"])
    params:
        input_dir=config["outdir"] + "/pmf_grn/experiments/full_prior/{dataset}/{c}",
        output_dir=config["outdir"] + "/pmf_grn/{dataset}/{c}"
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/pmf-grn/env"
    shell:
        """
        cd {config[algo][pmf_grn][dir_path]}
        python evaluate.py --exp-dirs {params.input_dir} \
        --output-dir {params.output_dir} --expression-names {wildcards.c} --num-sampling-iters 500
        """

####################
## Step 7: Get summary stats
####################
# rule pmf_grn_quad_get_summary_stat:
#     input:
#         A=expand("{outdir}/pmf_grn/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/inferred_consensus_grn.tsv", outdir=config["outdir"]),
#         P=expand("{outdir}/pmf_grn/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/inferred_consensus_tfa.tsv", outdir=config["outdir"]),
#         prior=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
#         true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
#         g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
#     output:
#         expand("{outdir}/pmf_grn/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
#         expand("{outdir}/pmf_grn/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
#     params:
#         outdir=expand("{outdir}/pmf_grn/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}", outdir=config["outdir"])
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
#         --transpose_inferred_P \
#         --transpose_inferred_A \
#         --inferred_A_ignore_first_row_col \
#         --inferred_P_ignore_first_row_col \
#         --not_scaling_inferred_A # this is important because the model includes sequencing depth as hidden var
#         """

# rule pmf_grn_get_clustering_accurarcy:
#     input: 
#         P=expand("{outdir}/pmf_grn/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/inferred_consensus_tfa.tsv", outdir=config["outdir"]),
#         gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
#     output:
#         outstat=expand("{outdir}/pmf_grn/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
#     conda:
#         "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/Seurat"
#     shell:
#         """
#         Rscript {config[scriptdir]}/pmf-grn/Seurat_clustering.R \
#         {input.P} \
#         {input.gt_celltype} \
#         {config[tflist]} \
#         $'\\t' \
#         {output.outstat}
#         """