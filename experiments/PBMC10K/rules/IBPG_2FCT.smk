rule run_ibpg_2fct_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/all_regulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/all_exps/all_exps.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/all_exps/all_priors.txt", datadir=config["datadir"]),
        t=expand("{datadir}/{{dataset}}/all_exps/celltype_tree_ancestor.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/ibpg_2fct/{{dataset}}/all_cell_results/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=cell_clusters),
        out_P=expand("{outdir}/ibpg_2fct/{{dataset}}/all_cell_results/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=cell_clusters),
        out_rel=expand("{outdir}/ibpg_2fct/{{dataset}}/all_cell_results/rep_0/relative_error.txt", outdir=config["outdir"]),
        outdir=directory(expand("{outdir}/ibpg_2fct/{{dataset}}/all_cell_results", outdir=config["outdir"]))
    params:
        exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
        num_rep=num_random_init
    threads: 20
    benchmark: 
        "benchmarks/ibpg_2fct/{dataset}.benchmark.txt"
    shell:
        """
        {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
        --max_iter {config[algo][ibpg_2fct][max_iter]} --init_m {config[algo][ibpg_2fct][init_m]} \
        --num_rep {params.num_rep} -t ibpg_2fct -l {config[algo][ibpg_2fct][l]} --lambda_p {config[algo][ibpg_2fct][lambda_p]} \
        --tree {input.t} --lambda_Cminus {config[algo][ibpg_2fct][lambda_Cminus]} --lambda_E {config[algo][ibpg_2fct][lambda_E]}
        """

# rule BH_run_ibpg_2fct_TFA_inference:
#     input:
#         r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/allregulators.txt", datadir=config["datadir"]),
#         d=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/all_exps/all_exps.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/all_exps/all_genes.txt", datadir=config["datadir"]),
#         p=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/all_exps/all_priors.txt", datadir=config["datadir"]),
#         t=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/all_exps/celltype_tree_ancestor.txt", datadir=config["datadir"]),
#     output:
#         out_A=expand("{outdir}/ibpg_2fct/processed_Buenrostro_Hematopoiesis/all_cell_results/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=BH_clusters),
#         out_P=expand("{outdir}/ibpg_2fct/processed_Buenrostro_Hematopoiesis/all_cell_results/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=BH_clusters),
#         out_rel=expand("{outdir}/ibpg_2fct/processed_Buenrostro_Hematopoiesis/all_cell_results/rep_0/relative_error.txt", outdir=config["outdir"]),
#         outdir=directory(expand("{outdir}/ibpg_2fct/processed_Buenrostro_Hematopoiesis/all_cell_results", outdir=config["outdir"]))
#     params:
#         exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
#         num_rep=num_random_init
#     threads: 20
#     benchmark: 
#         "benchmarks/ibpg_2fct/processed_Buenrostro_Hematopoiesis.benchmark.txt"
#     shell:
#         """
#         {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
#         --max_iter {config[algo][ibpg_2fct][max_iter]} --init_m {config[algo][ibpg_2fct][init_m]} \
#         --num_rep {params.num_rep} -t ibpg_2fct -l {config[algo][ibpg_2fct][l]} --lambda_p {config[algo][ibpg_2fct][lambda_p]} \
#         --tree {input.t} --lambda_Cminus {config[algo][ibpg_2fct][lambda_Cminus]} --lambda_E {config[algo][ibpg_2fct][lambda_E]}
#         """

# rule RFH_run_ibpg_2fct_TFA_inference:
#     input:
#         r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/allregulators.txt", datadir=config["datadir"]),
#         d=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/all_exps/all_exps.txt", datadir=config["datadir"]),
#         g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/all_exps/all_genes.txt", datadir=config["datadir"]),
#         p=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/all_exps/all_priors.txt", datadir=config["datadir"]),
#         t=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/all_exps/celltype_tree_ancestor.txt", datadir=config["datadir"]),
#     output:
#         out_A=expand("{outdir}/ibpg_2fct/processed_Ranzoni_fetal_hematopoiesis_fine_tree/all_cell_results/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=RFH_fine_clusters),
#         out_P=expand("{outdir}/ibpg_2fct/processed_Ranzoni_fetal_hematopoiesis_fine_tree/all_cell_results/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=RFH_fine_clusters),
#         out_rel=expand("{outdir}/ibpg_2fct/processed_Ranzoni_fetal_hematopoiesis_fine_tree/all_cell_results/rep_0/relative_error.txt", outdir=config["outdir"]),
#         outdir=directory(expand("{outdir}/ibpg_2fct/processed_Ranzoni_fetal_hematopoiesis_fine_tree/all_cell_results", outdir=config["outdir"]))
#     params:
#         exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
#         num_rep=num_random_init
#     threads: 20
#     benchmark: 
#         "benchmarks/ibpg_2fct/processed_Ranzoni_fetal_hematopoiesis_fine_tree.benchmark.txt"
#     shell:
#         """
#         {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
#         --max_iter {config[algo][ibpg_2fct][max_iter]} --init_m {config[algo][ibpg_2fct][init_m]} \
#         --num_rep {params.num_rep} -t ibpg_2fct -l {config[algo][ibpg_2fct][l]} --lambda_p {config[algo][ibpg_2fct][lambda_p]} \
#         --tree {input.t} --lambda_Cminus {config[algo][ibpg_2fct][lambda_Cminus]} --lambda_E {config[algo][ibpg_2fct][lambda_E]}
#         """

# as all A matrices from the previous rule follows the gene order of ibpg_2fct_g, we make filtered and reordered graphs based on the cell type specific gene order
# for this dataset this step should do nothing.
rule ibpg_2fct_filter_and_reorder_by_cell_type:
    input:
        A=expand("{outdir}/ibpg_2fct/{{dataset}}/all_cell_results/rep_0/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ibpg_2fct/{{dataset}}/all_cell_results/rep_0/{{c}}/P.txt", outdir=config["outdir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]), # as genes are the same for all cell types
        ibpg_2fct_g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/ibpg_2fct/{{dataset}}/{{c}}/rep_0/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/ibpg_2fct/{{dataset}}/{{c}}/rep_0/P.txt", outdir=config["outdir"]),
        out_rel=expand("{outdir}/ibpg_2fct/{{dataset}}/{{c}}/rep_0/relative_error.txt", outdir=config["outdir"])
    params:
        indir=expand("{outdir}/ibpg_2fct/{{dataset}}/all_cell_results", outdir=config["outdir"]),
        num_rep=num_random_init,
        outdir=expand("{outdir}/ibpg_2fct/{{dataset}}/{{c}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/IBPG_2FCT/filter_and_reorder_output_by_cell_type_genes.py --indir {params.indir} --outdir {params.outdir} \
            --in_gene_order {input.ibpg_2fct_g} --cell_type {wildcards.c} --cell_type_gene_order {input.g} --num_rep {params.num_rep}
        """


rule run_ibpg_2fct_500_1000_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/all_regulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/all_exps/all_exps.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/all_exps/all_priors.txt", datadir=config["datadir"]),
        t=expand("{datadir}/{{dataset}}/all_exps/celltype_tree_ancestor.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/all_cell_results/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=cell_clusters),
        out_P=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/all_cell_results/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=cell_clusters),
        out_rel=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/all_cell_results/rep_0/relative_error.txt", outdir=config["outdir"]),
        outdir=directory(expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/all_cell_results", outdir=config["outdir"]))
    params:
        exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
        num_rep=num_random_init
    threads: 20
    benchmark: 
        "benchmarks/ibpg_2fct_500_1000/{dataset}.benchmark.txt"
    shell:
        """
        {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
        --max_iter {config[algo][ibpg_2fct][max_iter]} --init_m {config[algo][ibpg_2fct][init_m]} \
        --num_rep {params.num_rep} -t ibpg_2fct -l 500 --lambda_p {config[algo][ibpg_2fct][lambda_p]} \
        --tree {input.t} --lambda_Cminus {config[algo][ibpg_2fct][lambda_Cminus]} --lambda_E 1000
        """


# as all A matrices from the previous rule follows the gene order of ibpg_2fct_g, we make filtered and reordered graphs based on the cell type specific gene order
# for this dataset this step should do nothing.
rule ibpg_2fct_500_1000_filter_and_reorder_by_cell_type:
    input:
        A=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/all_cell_results/rep_0/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/all_cell_results/rep_0/{{c}}/P.txt", outdir=config["outdir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]), # as genes are the same for all cell types
        ibpg_2fct_g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/{{c}}/rep_0/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/{{c}}/rep_0/P.txt", outdir=config["outdir"]),
        out_rel=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/{{c}}/rep_0/relative_error.txt", outdir=config["outdir"])
    params:
        indir=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/all_cell_results", outdir=config["outdir"]),
        num_rep=num_random_init,
        outdir=expand("{outdir}/ibpg_2fct_500_1000/{{dataset}}/{{c}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/IBPG_2FCT/filter_and_reorder_output_by_cell_type_genes.py --indir {params.indir} --outdir {params.outdir} \
            --in_gene_order {input.ibpg_2fct_g} --cell_type {wildcards.c} --cell_type_gene_order {input.g} --num_rep {params.num_rep}
        """
        

rule run_ibpg_2fct_500_0_TFA_inference:
    input:
        r=expand("{datadir}/{{dataset}}/all_regulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/all_exps/all_exps.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/all_exps/all_priors.txt", datadir=config["datadir"]),
        t=expand("{datadir}/{{dataset}}/all_exps/celltype_tree_ancestor.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/all_cell_results/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=cell_clusters),
        out_P=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/all_cell_results/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=cell_clusters),
        out_rel=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/all_cell_results/rep_0/relative_error.txt", outdir=config["outdir"]),
        outdir=directory(expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/all_cell_results", outdir=config["outdir"]))
    params:
        exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
        num_rep=num_random_init
    threads: 20
    benchmark: 
        "benchmarks/ibpg_2fct_500_0/{dataset}.benchmark.txt"
    shell:
        """
        {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
        --max_iter {config[algo][ibpg_2fct][max_iter]} --init_m {config[algo][ibpg_2fct][init_m]} \
        --num_rep {params.num_rep} -t ibpg_2fct -l 500 --lambda_p {config[algo][ibpg_2fct][lambda_p]} \
        --tree {input.t} --lambda_Cminus {config[algo][ibpg_2fct][lambda_Cminus]} --lambda_E 0
        """


# as all A matrices from the previous rule follows the gene order of ibpg_2fct_g, we make filtered and reordered graphs based on the cell type specific gene order
# for this dataset this step should do nothing.
rule ibpg_2fct_500_0_filter_and_reorder_by_cell_type:
    input:
        A=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/all_cell_results/rep_0/{{c}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/all_cell_results/rep_0/{{c}}/P.txt", outdir=config["outdir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]), # as genes are the same for all cell types
        ibpg_2fct_g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/{{c}}/rep_0/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/{{c}}/rep_0/P.txt", outdir=config["outdir"]),
        out_rel=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/{{c}}/rep_0/relative_error.txt", outdir=config["outdir"])
    params:
        indir=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/all_cell_results", outdir=config["outdir"]),
        num_rep=num_random_init,
        outdir=expand("{outdir}/ibpg_2fct_500_0/{{dataset}}/{{c}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/IBPG_2FCT/filter_and_reorder_output_by_cell_type_genes.py --indir {params.indir} --outdir {params.outdir} \
            --in_gene_order {input.ibpg_2fct_g} --cell_type {wildcards.c} --cell_type_gene_order {input.g} --num_rep {params.num_rep}
        """