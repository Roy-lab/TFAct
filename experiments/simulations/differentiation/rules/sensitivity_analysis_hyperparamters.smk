


rule ibpg_2fct_hyper:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/{{TFAsimrep}}/{{exp}}_all_exps.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_genes.txt", datadir=config['datadir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_priors.txt", datadir=config["datadir"]),
        t=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/cell_lineage.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=prefixes),
        out_P=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=prefixes),
        out_rel=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/relative_error.txt", outdir=config["outdir"]),
        outdir=directory(expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"]))
    params:
        exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
        num_rep=num_random_init,
        lambda_p = lambda_p,
        lambda_Cminus = lambda_Cminus
    threads: 10
    benchmark: 
        "benchmarks/ibpg_2fct/{l1}_{l_tree}/{edgeprob}/{ldak}/{TFAsimrep}/{exp}.benchmark.txt"
    shell:
        """
        {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
        --max_iter 500 --init_m {config[algo][ibpg_2fct][init_m]} \
        --num_rep {params.num_rep} -t ibpg_2fct -l {wildcards.l1} --lambda_p {params.lambda_p} \
        --tree {input.t} --lambda_Cminus {params.lambda_Cminus} --lambda_E {wildcards.l_tree}
        """

# as all A matrices from the previous rule follows the gene order of ibpg_2fct_g, we make filtered and reordered graphs based on the cell type gene order
# this step is currently inefficient. Can reduce the number of jobs by merging all cell types into one job. 
rule ibpg_2fct_reorder_hyper:
    input:
        A=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{{pre}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{{pre}}/P.txt", outdir=config["outdir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        ibpg_2fct_g=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_genes.txt", datadir=config['datadir'])
    output:
        out_A=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/P.txt", outdir=config["outdir"]),
        out_rel=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/relative_error.txt", outdir=config["outdir"])
    params:
        indir=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"]),
        num_rep=num_random_init,
        outdir=expand("{outdir}/ct_algos/ibpg_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/IBPG_2FCT/filter_and_reorder_output_by_cell_type_genes.py --indir {params.indir} --outdir {params.outdir} \
            --in_gene_order {input.ibpg_2fct_g} --cell_type {wildcards.pre} --cell_type_gene_order {input.g} --num_rep {params.num_rep}
        """


rule kl_nmf_2fct_hyper:
    input:
        r = config["tflist"],
        d=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/{{TFAsimrep}}/{{exp}}_all_exps.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_genes.txt", datadir=config['datadir']),
        p=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_priors.txt", datadir=config["datadir"]),
        t=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/cell_lineage.txt", datadir=config["datadir"])
    output:
        out_A=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{pre}/A.txt", outdir=config["outdir"], pre=prefixes),
        out_P=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{pre}/P.txt", outdir=config["outdir"], pre=prefixes),
        out_rel=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/relative_error.txt", outdir=config["outdir"]),
        outdir=directory(expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"]))
    params:
        exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/IBPG/tetraNMF",
        num_rep=num_random_init,
        lambda_p = lambda_p,
        lambda_Cminus = lambda_Cminus
    threads: 10
    benchmark: 
        "benchmarks/kl_nmf_2fct/{l1}_{l_tree}/{edgeprob}/{ldak}/{TFAsimrep}/{exp}.benchmark.txt"
    shell:
        """
        {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
        --max_iter 500 --init_m {config[algo][kl_nmf_2fct][init_m]} \
        --num_rep {params.num_rep} -t KL_NMF_2fct -l {wildcards.l1} --lambda_p {params.lambda_p} \
        --tree {input.t} --lambda_Cminus {params.lambda_Cminus} --lambda_E {wildcards.l_tree}
        """

# as all A matrices from the previous rule follows the gene order of kl_nmf_2fct_g, we make filtered and reordered graphs based on the cell type specific gene order
rule kl_NMF_2fct_reorder_hyper:
    input:
        A=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{{pre}}/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}/rep_0/{{pre}}/P.txt", outdir=config["outdir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
        kl_nmf_2fct_g=expand("{datadir}/{{edgeprob}}/{{ldak}}/all_exps/all_genes.txt", datadir=config['datadir'])
    output:
        out_A=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/P.txt", outdir=config["outdir"]),
        out_rel=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/relative_error.txt", outdir=config["outdir"])
    params:
        indir=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/all_cell_results/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"]),
        num_rep=num_random_init,
        outdir=expand("{outdir}/ct_algos/kl_nmf_2fct/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}", outdir=config["outdir"])
    shell:
        """
        python {config[scriptdir]}/IBPG_2FCT/filter_and_reorder_output_by_cell_type_genes.py --indir {params.indir} --outdir {params.outdir} \
            --in_gene_order {input.kl_nmf_2fct_g} --cell_type {wildcards.pre} --cell_type_gene_order {input.g} --num_rep {params.num_rep}
        """

####################
## Get summary stats
####################
rule ct_stability_selection_get_mean:
    input:
        rel=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/relative_error.txt", outdir=config["outdir"]),
        in_A=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/A.txt", outdir=config["outdir"]),
        in_P=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}/rep_0/P.txt", outdir=config["outdir"])
    output:
        out_A=expand("{outdir}/ct_algos/{{algo}}/{{l1,[\.0-9]+}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/A.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/ct_algos/{{algo}}/{{l1,[\.0-9]+}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/P.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(ct_algos) 
    params:
        indir_suffix="ct_algos/{algo}/{l1}_{l_tree}/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}",
        output_precision = 7,
        file_pattern_A = "rep_*/A.txt",
        file_pattern_P = "rep_*/P.txt",
        frac_ignored= 0.25
    shell:
        """
        python {config[scriptdir]}/get_stability_selection_mean.py --A_file_pattern '{config[outdir]}/{params.indir_suffix}/{params.file_pattern_A}' \
        --P_file_pattern '{config[outdir]}/{params.indir_suffix}/{params.file_pattern_P}' --frac_ignored {params.frac_ignored} \
        --output_precision {params.output_precision} --out_A {output.out_A} --out_P {output.out_P}
        """

rule ct_get_summary_stat:
    input:
        A=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(ct_algos) 
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # TODO change this
    params:
        outdir=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}", outdir=config["outdir"])
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
        --transpose_inferred_A
        """

rule ct_collect_F_scores:
    input:
        expand("{outdir}/ct_algos/{algo}/{l1}_{l_tree}/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}_stability/stats_{top_e}/stat.txt", outdir=config["outdir"], \
        algo=ct_algos, l1=l1_choices, l_tree=lineage_reg_choices, edgeprob=edge_adding_prob,ldak=lda_k,pre=prefixes,TFAsimrep=TFA_sim_rep, exp=exp_mat_w_dropout, top_e=top_edges)
    output:
        expand("{outdir}/analysis/ct/Fscores.txt", outdir=config["outdir"])
    params:
        stat_path=[config["outdir"]+ '/ct_algos/'+alg + '/*_*/*/*/*/*/*/stats_*/stat.txt' for alg in ct_algos]
    shell: # note: modify below to include more algorithms
        """
        awk 'BEGIN{{OFS="\\t"}}NR==1{{$1=$1; printf "file_name\\t%s\\n", $0}}FNR==2{{printf "%s\\t",  FILENAME}} FNR>1{{print}}' {params.stat_path} | \
        awk -F'\\t' -v regex_pattern='{config[outdir]}/ct_algos/([^/]+|baselines/[^/]+)/([\\.0-9]+)_([\\.0-9]+)/(0.[0-9]+)/(k[0-9]+)/([^/]+)/(rep[0-9]+)/([^/]+)_stability/stats_([0-9]+)'  'BEGIN{{OFS="\\t"}} NR==1{{printf "algo\\tl_prior\\tl_sim\\tedge_addition_prob\\tlda_topics\\tcell\\tTFA_rep\\texp_dropout\\tFscore_top_edge\\t"; 
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF }}NR>1{{
            $1=$1; 
            match($1, regex_pattern, arr);
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t", arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7], arr[8], arr[9];
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF 
        }}' > {output}
        """

rule ct_get_clustering_accurarcy:
    input: 
        P=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/P.txt", outdir=config["outdir"]),
        gt_celltype=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/cell_celltype_assignments.txt", datadir=config["datadir"]),
    output:
        outstat=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(ct_algos) 
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

rule ct_collect_clust_accurarcy:
    input:
        expand("{outdir}/ct_algos/{algo}/{l1}_{l_tree}/{edgeprob}/{ldak}/{pre}/{TFAsimrep}/{exp}_stability/downstream_analysis/clust_stats.txt", outdir=config["outdir"], \
        algo=ct_algos, l1=l1_choices, l_tree=lineage_reg_choices, edgeprob=edge_adding_prob,ldak=lda_k,pre=prefixes,TFAsimrep=TFA_sim_rep, exp=exp_mat_w_dropout, top_e=top_edges)
    output:
        expand("{outdir}/analysis/ct/clustering_stats.txt", outdir=config["outdir"])
    params:
        stat_path=[config["outdir"]+ '/ct_algos/'+alg + '/*_*/*/*/*/*/*/downstream_analysis/clust_stats.txt' for alg in ct_algos]
    shell: # note: modify below to include more algorithms
        """
        awk 'BEGIN{{OFS="\\t"}}NR==1{{$1=$1; printf "file_name\\t%s\\n", $0}}FNR==2{{printf "%s\\t",  FILENAME}} FNR>1{{print}}' {params.stat_path} | \
        awk -F'\\t' -v regex_pattern='{config[outdir]}/ct_algos/([^/]+|baselines/[^/]+)/([\\.0-9]+)_([\\.0-9]+)/(0.[0-9]+)/(k[0-9]+)/([^/]+)/(rep[0-9]+)/([^/]+)_stability/downstream_analysis'  'BEGIN{{OFS="\\t"}} NR==1{{printf "algo\\tl_prior\\tl_sim\\tedge_addition_prob\\tlda_topics\\tcell\\tTFA_rep\\texp_dropout\\t"; 
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF }}NR>1{{
            $1=$1; 
            match($1, regex_pattern, arr);
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t", arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7], arr[8];
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF 
        }}' > {output}
        """


### 
# AUPR and reg specific AUPR
###

rule ct_get_AUPR:
    input:
        A=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/P.txt", outdir=config["outdir"]),
        gs=expand("{genedir}/{{edgeprob}}/{{pre}}/nw.tsv", genedir=config["genedir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/AUPR_TF.txt", outdir=config["outdir"]),
    wildcard_constraints:
        algo="|".join(ct_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR/", outdir=config["outdir"]),
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

rule ct_get_AUPR_prior:
    input:
        A=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}_prior.txt", datadir=config["datadir"]),
        r=expand("{tflist}", tflist=config["tflist"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{genedir}/{{edgeprob}}/{{pre}}/genes.txt", genedir=config['genedir']),
    output:
        expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/AUPR_TF.txt", outdir=config["outdir"]),
    wildcard_constraints:
        algo="|".join(ct_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/ct_algos/{{algo}}/{{l1}}_{{l_tree}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/AUPR_prior/", outdir=config["outdir"]),
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