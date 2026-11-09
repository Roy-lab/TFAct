# rule kl_NMF_2fct_make_input:
#     input:
#         in_gene_dir=expand("{genedir}/{{edgeprob}}", genedir=config["genedir"]), 
#         in_exp_prior_dir=expand("{datadir}/{{edgeprob}}/{{ldak}}", datadir=config["datadir"]),
#     output:
#         d=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/{TFAsimrep}/{exp}_all_exps.txt", \
#             datadir=config["datadir"], TFAsimrep=TFA_sim_rep, exp=exp_mat_w_dropout),
#         g=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_genes.txt", datadir=config['datadir']),
#         p=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_priors.txt", datadir=config["datadir"]),
#     params:
#         outdir=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps", datadir=config["datadir"]),
#         TFA_sim_rep=lambda wildcards: ','.join(f"{_}" for _ in TFA_sim_rep),
#         exp_name=lambda wildcards: ','.join(f"{_}" for _ in exp_mat_w_dropout),
#         cell_types=lambda wildcards: ','.join(f"{_}" for _ in prefixes)
#     shell:
#         """
#         bash {config[scriptdir]}/IBPG_2FCT/make_inputs_v2.sh {input.in_gene_dir} {input.in_exp_prior_dir} {params.outdir} \
#             {params.TFA_sim_rep} {params.exp_name} {params.cell_types} {params.cell_types}
#         """

rule run_kl_NMF_2fct_TFA_inference:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/{{TFAsimrep}}/{{exp}}_all_exps.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_genes.txt", datadir=config['datadir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_priors.txt", datadir=config["datadir"]),
        t=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/cell_lineage.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=prefixes),
        out_P=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=prefixes),
        out_rel=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/relative_error.txt", outdir=config["outdir"]),
        outdir=directory(expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"]))
    params:
        exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
        num_rep=num_random_init
    threads: 10
    benchmark: 
        "benchmarks/kl_nmf_2fct/{edgeprob}/{ldak}/{TFAsimrep}/{exp}.benchmark.txt"
    shell:
        """
        {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
        --max_iter {config[algo][kl_nmf_2fct][max_iter]} --init_m {config[algo][kl_nmf_2fct][init_m]} \
        --num_rep {params.num_rep} -t KL_NMF_2fct -l {config[algo][kl_nmf_2fct][l]} --lambda_p {config[algo][kl_nmf_2fct][lambda_p]} \
        --tree {input.t} --lambda_Cminus {config[algo][kl_nmf_2fct][lambda_Cminus]} --lambda_E {config[algo][kl_nmf_2fct][lambda_E]}
        """

# as all A matrices from the previous rule follows the gene order of kl_nmf_2fct_g, we make filtered and reordered graphs based on the cell type specific gene order
rule kl_NMF_2fct_filter_and_reorder_by_cell_type:
    input:
        A=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{{pre}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{{pre}}/P.txt", outdir=config["outdir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        kl_nmf_2fct_g=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_genes.txt", datadir=config['datadir'])
    output:
        out_A=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/P.txt", outdir=config["outdir"]),
        out_rel=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/relative_error.txt", outdir=config["outdir"])
    params:
        indir=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"]),
        num_rep=num_random_init,
        outdir=expand("{outdir}/kl_nmf_2fct/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/IBPG_2FCT/filter_and_reorder_output_by_cell_type_genes.py --indir {params.indir} --outdir {params.outdir} \
            --in_gene_order {input.kl_nmf_2fct_g} --cell_type {wildcards.pre} --cell_type_gene_order {input.g} --num_rep {params.num_rep}
        """
        