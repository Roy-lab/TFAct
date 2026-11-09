rule FBS_algo_get_summary_stat:
    input:
        A=expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/{{algo}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(stability_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
        tf_by_gene=lambda wildcards: ' '.join(f"'{_}'" for _ in tf_by_gene_A_algos),
        gene_by_tf=lambda wildcards: ' '.join(f"'{_}'" for _ in gene_by_tf_A_algos)
    shell: # TODO add true_P;
        """
        tf_by_gene=( {params.tf_by_gene} )
        gene_by_tf=( {params.gene_by_tf} )
        if [[ " ${{gene_by_tf[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --transpose_inferred_A
        elif [[ " ${{tf_by_gene[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e}
        else # edge list
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        fi
        """


rule FBS_collect_F_scores:
    input:
        expand("{outdir}/{algo}/processed_FBS_reprogramming/{c}_stability/{gs}/stats_{top_e}/stat.txt", outdir=config["outdir"], \
        algo=algos, c=FBS_clusters, gs=mouse_gold_standards, top_e=top_edges)
    output:
        expand("{outdir}/processed_FBS_reprogramming/analysis/Fscores.txt", outdir=config["outdir"])
    params:
        stat_path=[config["outdir"]+ '/'+alg + '/processed_FBS_reprogramming/*/*/stats_*/stat.txt' for alg in algos]
    shell: # note: modify below to include more algorithms
        """
        awk 'BEGIN{{OFS="\\t"}}NR==1{{$1=$1; printf "file_name\\t%s\\n", $0}}FNR==2{{printf "%s\\t",  FILENAME}} FNR>1{{print}}' {params.stat_path} | \
        awk -F'\\t' -v regex_pattern='{config[outdir]}/([^/]+|baselines/[^/]+)/processed_FBS_reprogramming/([^/]+)_stability/([^/]+)/stats_([0-9]+)'  'BEGIN{{OFS="\\t"}} NR==1{{printf "algo\\tcell\\tgold_standard\\tFscore_top_edge\\t"; 
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF }}NR>1{{
            $1=$1; 
            match($1, regex_pattern, arr);
            printf "%s\\t%s\\t%s\\t%s\\t", arr[1], arr[2], arr[3], arr[4];
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF 
        }}' > {output}
        """


rule BH_algo_get_summary_stat:
    input:
        A=expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/{{algo}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(stability_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
        tf_by_gene=lambda wildcards: ' '.join(f"'{_}'" for _ in tf_by_gene_A_algos),
        gene_by_tf=lambda wildcards: ' '.join(f"'{_}'" for _ in gene_by_tf_A_algos)
    shell: # TODO add true_P;
        """
        tf_by_gene=( {params.tf_by_gene} )
        gene_by_tf=( {params.gene_by_tf} )
        if [[ " ${{gene_by_tf[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --transpose_inferred_A
        elif [[ " ${{tf_by_gene[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e}
        else # edge list
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        fi
        """


rule BH_collect_F_scores:
    input:
        expand("{outdir}/{algo}/processed_Buenrostro_Hematopoiesis/{c}_stability/{gs}/stats_{top_e}/stat.txt", outdir=config["outdir"], \
        algo=algos, c=BH_clusters, gs=BH_gold_standards, top_e=top_edges)
    output:
        expand("{outdir}/processed_Buenrostro_Hematopoiesis/analysis/Fscores.txt", outdir=config["outdir"])
    params:
        stat_path=[config["outdir"]+ '/'+alg + '/processed_Buenrostro_Hematopoiesis/*/*/stats_*/stat.txt' for alg in algos]
    shell: # note: modify below to include more algorithms
        """
        awk 'BEGIN{{OFS="\\t"}}NR==1{{$1=$1; printf "file_name\\t%s\\n", $0}}FNR==2{{printf "%s\\t",  FILENAME}} FNR>1{{print}}' {params.stat_path} | \
        awk -F'\\t' -v regex_pattern='{config[outdir]}/([^/]+|baselines/[^/]+)/processed_Buenrostro_Hematopoiesis/([^/]+)_stability/([^/]+)/stats_([0-9]+)'  'BEGIN{{OFS="\\t"}} NR==1{{printf "algo\\tcell\\tgold_standard\\tFscore_top_edge\\t"; 
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF }}NR>1{{
            $1=$1;
            match($1, regex_pattern, arr);
            printf "%s\\t%s\\t%s\\t%s\\t", arr[1], arr[2], arr[3], arr[4];
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF 
        }}' > {output}
        """

rule RFH_algo_get_summary_stat:
    input:
        A=expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/P.txt", outdir=config["outdir"]),
        prior=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_network.txt", datadir=config["datadir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # expand("{outdir}/{{algo}}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/{{exp}}_stability/stats_{{top_e}}/tf_specifics.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}/stat.txt", outdir=config["outdir"])
    wildcard_constraints:
        algo="|".join(stability_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
        tf_by_gene=lambda wildcards: ' '.join(f"'{_}'" for _ in tf_by_gene_A_algos),
        gene_by_tf=lambda wildcards: ' '.join(f"'{_}'" for _ in gene_by_tf_A_algos)
    shell: # TODO add true_P;
        """
        tf_by_gene=( {params.tf_by_gene} )
        gene_by_tf=( {params.gene_by_tf} )
        if [[ " ${{gene_by_tf[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --transpose_inferred_A
        elif [[ " ${{tf_by_gene[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e}
        else # edge list
            python {config[scriptdir]}/compare_networks_and_TFA_v2.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --prior_A {input.prior} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        fi
        """

# both BH and RFH are from human samples so use the same gold standard
rule RFH_collect_F_scores:
    input:
        expand("{outdir}/{algo}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{c}_stability/{gs}/stats_{top_e}/stat.txt", outdir=config["outdir"], \
        algo=algos, c=RFH_fine_clusters, gs=BH_gold_standards, top_e=top_edges)
    output:
        expand("{outdir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/analysis/Fscores.txt", outdir=config["outdir"])
    params:
        stat_path=[config["outdir"]+ '/'+alg + '/processed_Ranzoni_fetal_hematopoiesis_fine_tree/*/*/stats_*/stat.txt' for alg in algos]
    shell: # note: modify below to include more algorithms
        """
        awk 'BEGIN{{OFS="\\t"}}NR==1{{$1=$1; printf "file_name\\t%s\\n", $0}}FNR==2{{printf "%s\\t",  FILENAME}} FNR>1{{print}}' {params.stat_path} | \
        awk -F'\\t' -v regex_pattern='{config[outdir]}/([^/]+|baselines/[^/]+)/processed_Ranzoni_fetal_hematopoiesis_fine_tree/([^/]+)_stability/([^/]+)/stats_([0-9]+)'  'BEGIN{{OFS="\\t"}} NR==1{{printf "algo\\tcell\\tgold_standard\\tFscore_top_edge\\t"; 
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF }}NR>1{{
            $1=$1;
            match($1, regex_pattern, arr);
            printf "%s\\t%s\\t%s\\t%s\\t", arr[1], arr[2], arr[3], arr[4];
            for (i=2; i<NF; i++) printf $i "\\t"; print $NF 
        }}' > {output}
        """

###
### AUPR and predictable TFs
###


rule FBS_algo_get_AUPR:
    input:
        A=expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/mouse_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_FBS_reprogramming/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        # AUPR has nothing to do with top edges. this is just to avoid snakemake bugs
        expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    wildcard_constraints:
        algo="|".join(stability_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
        tf_by_gene=lambda wildcards: ' '.join(f"'{_}'" for _ in tf_by_gene_A_algos),
        gene_by_tf=lambda wildcards: ' '.join(f"'{_}'" for _ in gene_by_tf_A_algos)
    shell: # TODO add true_P;
        """
        tf_by_gene=( {params.tf_by_gene} )
        gene_by_tf=( {params.gene_by_tf} )
        if [[ " ${{gene_by_tf[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --transpose_inferred_A
        elif [[ " ${{tf_by_gene[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e}
        else # edge list
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        fi
        """



rule BH_algo_get_AUPR:
    input:
        A=expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Buenrostro_Hematopoiesis/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    wildcard_constraints:
        algo="|".join(stability_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/{{algo}}/processed_Buenrostro_Hematopoiesis/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
        tf_by_gene=lambda wildcards: ' '.join(f"'{_}'" for _ in tf_by_gene_A_algos),
        gene_by_tf=lambda wildcards: ' '.join(f"'{_}'" for _ in gene_by_tf_A_algos)
    shell: # TODO add true_P;
        """
        tf_by_gene=( {params.tf_by_gene} )
        gene_by_tf=( {params.gene_by_tf} )
        if [[ " ${{gene_by_tf[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --transpose_inferred_A
        elif [[ " ${{tf_by_gene[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e}
        else # edge list
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        fi
        """

rule RFH_algo_get_AUPR:
    input:
        A=expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/A.txt", outdir=config["outdir"]),
        P=expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/P.txt", outdir=config["outdir"]),
        gs=expand("{datadir}/Goldstandards/human_goldstandard/{{gs}}.txt", datadir=config["datadir"]),
        r=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allregulators.txt", datadir=config["datadir"]),
        # true_P=expand("{datadir}/{{edgeprob}}/{{ldak}}/{{pre}}/{{TFAsimrep}}/P.txt", datadir=config["datadir"]),
        g=expand("{datadir}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/inputdata/{{c}}_allGenes.txt", datadir=config["datadir"]),
    output:
        expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}/A_inferred_top_{{top_e}}.txt", outdir=config["outdir"]),
        expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
    wildcard_constraints:
        algo="|".join(stability_algos)
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/scanpy" # change this
    params:
        outdir=expand("{outdir}/{{algo}}/processed_Ranzoni_fetal_hematopoiesis_fine_tree/{{c}}_stability/{{gs}}/stats_{{top_e}}", outdir=config["outdir"]),
        tf_by_gene=lambda wildcards: ' '.join(f"'{_}'" for _ in tf_by_gene_A_algos),
        gene_by_tf=lambda wildcards: ' '.join(f"'{_}'" for _ in gene_by_tf_A_algos)
    shell: # TODO add true_P;
        """
        tf_by_gene=( {params.tf_by_gene} )
        gene_by_tf=( {params.gene_by_tf} )
        if [[ " ${{gene_by_tf[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --transpose_inferred_A
        elif [[ " ${{tf_by_gene[*]}} " =~ [[:space:]]{wildcards.algo}[[:space:]] ]]; then
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e}
        else # edge list
            python {config[scriptdir]}/compare_networks_and_TFA_AUPR.py \
            --outdir {params.outdir} \
            --inferred_A {input.A}  \
            --inferred_P {input.P} \
            --true_A {input.gs} \
            --tf_lst {input.r} \
            --gene {input.g} \
            --top {wildcards.top_e} \
            --A_edge_lst
        fi
        """



# rule FBS_algo_get_predictable_TF:
#     input:
#         fl=expand("{outdir}/{{algo}}/processed_FBS_reprogramming/{{c}}_stability/{{gs}}/stats_{{top_e}}/Enr_filelist_top_{{top_e}}.txt", outdir=config["outdir"]),
#     output:

#     params: # TODO
#         gsl_lib="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/brand_multiome_project/scMTNI_with_TFA/gsl/lib",
#         exe="/mnt/dv/wid/projects2/Roy-regnet-inference/yeast_regnet/yeast_regnet_from_scratch/may_6/test_new_ptf/comparetargetsets_allnets/compareTgtSetAcrossDataset"
    