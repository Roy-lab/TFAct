###
# Step 1: make scMTNI inputs
###

# should have be more space efficient. The following steps follow the scMTNI convention.
rule make_scMTNI_w_TFA_prior:
    input:
        p=expand("{datadir}/{{dataset}}/inputdata/{{c}}_network.txt",
                 datadir=config["datadir"]),
    output:
        p=expand("{datadir}/{{dataset}}/inputdata/scMTNI/prior/{{c}}_prior.txt",
                 datadir=config["datadir"]),
    shell:
        """
        # Add _cluster back to regulator and target. Make TFA prior as well
        awk -F'\\t' -v C="_{wildcards.c}" -v OFS='\\t' \
            '{{$1=$1 C; $2=$2 C; print; sub(C "$", "_TFA" C, $1); print}}' \
            {input.p} > {output.p}
        """


rule make_scMTNI_w_TFA_concatenated_expression:
    input:
        d=expand("{datadir}/{{dataset}}/inputdata/{{c}}.table",
                 datadir=config["datadir"]),
        c=expand("{datadir}/{{dataset}}/inputdata/{{c}}_cell_names.txt",
                 datadir=config["datadir"]),
        tf_lst=expand("{datadir}/{{dataset}}/inputdata/allregulators.txt", \
            datadir=config["datadir"], subsamp=str(num_subsamples-1)),
        tfa=expand("{TFAdir}/{{algo}}/{{dataset}}/{{c}}_stability/P.txt", \
            TFAdir=config["TFAdir"]),
    output:
        d=expand("{datadir}/{{dataset}}/inputdata/scMTNI/data/{{algo}}/{{c}}.table",
                 datadir=config["datadir"]),
    shell:
        """
        {{
            printf 'Gene\\t'
            paste -sd '\\t' {input.c}

            awk -F'\\t' -v C="_{wildcards.c}" -v OFS='\\t' \
                '{{$1=$1 C; print}}' {input.d}

            paste \
                <(awk -v C="_TFA_{wildcards.c}" '{{print $0 C}}' {input.tf_lst}) \
                <(awk -v OFS='\\t' '{{$1=$1; print}}' {input.tfa})
        }} > {output.d}
        """



rule extract_union_gene_set:
    input:
        d=lambda wildcards: expand(
            "{datadir}/{dataset}/inputdata/scMTNI/data/{algo}/{c}.table",
            datadir=config["datadir"],
            dataset=wildcards.dataset,
            algo=wildcards.algo,
            c=ds_cluster_map[wildcards.dataset],
        )
    output:
        g=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/Union.txt",
    params:
        outdir=lambda wildcards: (
            f"{config['datadir']}/{wildcards.dataset}/"
            f"inputdata/scMTNI/data/{wildcards.algo}"
        ),
        celltype_order=lambda wildcards: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[wildcards.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
    shell:
        """
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/Extracting_Union_geneset.sh \
        {params.outdir} \
        {params.celltype_order}
        """

rule make_OGID:
    input:
        g=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/Union.txt",
    output:
        ogid=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/OGID.txt",
    params:
        outdir=lambda wildcards: (
            f"{config['datadir']}/{wildcards.dataset}/"
            f"inputdata/scMTNI/data/{wildcards.algo}"
        ),
        celltype_order=lambda wildcards: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[wildcards.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
    shell:
        """
        cd {params.outdir} && \
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/OGID_making.sh \
        {params.celltype_order}
        """

rule make_allregulators_TFA_combined:
    input:
        r=lambda wildcards: expand(
            "{datadir}/{dataset}/inputdata/{c}_allregulators.txt",
            datadir=config["datadir"],
            dataset=celltype_order_dataset_map[wildcards.dataset],
            c=ds_cluster_map[wildcards.dataset],
        )
    output:
        r=config["datadir"] + "/{dataset}/inputdata/scMTNI/tfs/allregulators_TFA_combined.txt",
    params:
        clusters=lambda wildcards: "|".join(ds_cluster_map[wildcards.dataset]),
    shell:
        """
        mkdir -p $(dirname {output.r})

        cat {input.r} \
        | sed -E 's/_({params.clusters})$//' \
        | sort -u \
        | awk '{{print; print $0 "_TFA"}}' \
        > {output.r}
        """

rule create_scMTNI_file_lists:
    input:
        g=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/Union.txt",
        r=config["datadir"] + "/{dataset}/inputdata/scMTNI/tfs/allregulators_TFA_combined.txt",
    output:
        done=touch(config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/filelists/filelist.txt"),
    params:
        datadir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
        ),
        celltype_order=lambda w: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[w.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
    shell:
        """
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/create_file_lists.sh \
            {params.datadir} \
            {params.celltype_order}
        """

# note: change the module load command to an environment for the corresponding scripts in scMTNI under the conda directive.
rule run_scMTNI_other_input_file_creation:
    input:
        filelists=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/filelists/filelist.txt",
        tfs=config["datadir"] + "/{dataset}/inputdata/scMTNI/tfs/allregulators_TFA_combined.txt",
    output:
        done=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/testdata_ogids.txt",
    params:
        datadir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
        ),
        outdir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/results/{w.algo}"
        ),
        filelists=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}/filelists"
        ),
    threads: 5
    shell:
        """
        mkdir -p {params.outdir}

        module load bioconda-2024.02 && \
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/scMTNI.sh \
            {params.datadir} \
            {params.outdir} \
            {params.filelists} \
            {input.tfs} \
            {config[scriptdir]}/downstream_network_inference/scMTNI/PreparescMTNIinputfiles.py

        module unload bioconda-2024.02
        """

# create 50 subsamples
rule subsample_scMTNI_data:
    input:
        filelist=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/filelists/filelist.txt",
    output:
        done=touch(config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/subsample/done.txt"),
    params:
        indir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
        ),
        nseed=scMTNI_num_stability_subsamples,
    threads: 5
    shell:
        """
        python {config[scriptdir]}/downstream_network_inference/scMTNI/Datasubsample_sc_merged.py \
            --filelist {input.filelist} \
            --indir {params.indir} \
            --nseed {params.nseed}
        """

rule tar_scMTNI_prior:
    input:
        p=lambda w: expand(
            "{datadir}/{dataset}/inputdata/scMTNI/prior/{c}_prior.txt",
            datadir=config["datadir"],
            dataset=w.dataset,
            c=ds_cluster_map[w.dataset],
        )
    output:
        tar=config["datadir"] + "/{dataset}/inputdata/scMTNI/prior.tar.gz",
    params:
        indir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI"
        ),
    shell:
        """
        tar -czf {output.tar} -C {params.indir} prior
        """

rule make_scMTNI_condor_config:
    input:
        celltype_order=lambda w: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[w.dataset]}/"
            f"inputdata/celltype_order.txt"
        )
    output:
        configfile=config["datadir"] + "/{dataset}/inputdata/scMTNI/data/{algo}/condor_config.txt",
    shell:
        r"""
        awk -v OFS='\t' '
        NF {{
            c=$1
            print c, c, "results/" c, c "_allregulators.txt", c "_allGenes.txt", "prior/" c "_prior.txt"
        }}' {input.celltype_order} > {output.configfile}
        """


###
# Step 2: after running scMTNI, get the outputs and postprocess using kmeans and lda
###

rule organize_condor_outputs:
    output:
        expand("{outdir}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/randseed0/ogid0/results.tar.gz", outdir=config['outdir'])
    params:
        s=scMTNI_num_stability_subsamples,
        result_dir=expand("{outdir}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/", outdir=config["outdir"]),
        ogid_dir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}/ogids/"
        ),
    threads: 10
    shell:
        """
        num_ogids=$(ls {params.ogid_dir} | wc -l)
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/organize_condor_outputs_scMTNI.sh {params.result_dir} {params.result_dir} {params.s} $num_ogids
        """

# consider tarring up the randseed after this step to reduce inodes
rule Makeconsensusnetwork: # TODO: uncomment input once done
    # input:
    #     expand("{outdir}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/randseed0/ogid0/results.tar.gz", outdir=config['outdir'])
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.consensus.done")
    params:
        s=scMTNI_num_stability_subsamples,
        maxR=maxReg,
        ogid_dir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}/ogids/"
        ),
        indir=expand("{outdir}/algos/scMTNI/{{algo}}/{{dataset}}/4_1", outdir = config['outdir']),
        estEdge=estEdge,
        celltype_order=lambda wildcards: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[wildcards.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
    threads: 10
    shell:
        """
        num_ogids=$(ls {params.ogid_dir} | wc -l)
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/Makeconsensusnetwork.sh {params.s} {params.maxR} {params.celltype_order} \
        {params.indir} $num_ogids 0 {config[scriptdir]}/downstream_network_inference/scMTNI/Makeconsensusnetwork1.sh {params.estEdge}
        """

rule netfilterlowexpression:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.consensus.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.netfilter.done")
    params:
        celltype_order=lambda w: f"{config['datadir']}/{celltype_order_dataset_map[w.dataset]}/inputdata/celltype_order.txt",
        cvnpath=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/",
        indir=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}/"
    log:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/netfilterlowexpression.log"
    shell:
        """
        module load bioconda-2024.02
        Rscript --vanilla {config[scriptdir]}/downstream_network_inference/scMTNI/netfilterlowexpression.R \
            {params.celltype_order} {params.cvnpath} {params.indir} 1 > {log} 2>&1
        """

##it will filter _TFA targets (and _TFA regulators if they are expressed in less than 5 cells in the cell cluster)
rule extractTopconsensusedges_top:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.netfilter.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.top_{{top}}.done")
    log:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/logs_extractTopconsensusedges/top_{{top}}.log"
    params:
        analysis=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/",
        cell_order=lambda w: f"{config['datadir']}/{celltype_order_dataset_map[w.dataset]}/inputdata/celltype_order.txt",
        datadir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
        ),
    shell:
        """
        module load bioconda-2024.02
        Rscript --vanilla {config[scriptdir]}/downstream_network_inference/scMTNI/extractTopconsensusedges_scMTNI.R \
            {params.analysis} top {wildcards.top} {params.cell_order} {params.datadir} > {log} 2>&1
        """

##
# Kmeans
##
##it will filter _TFA targets (and _TFA regulators if they are expressed in less than 5 cells in the cell cluster)
rule extractTopconsensusedges_thr:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.netfilter.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.cf_{{thr}}.done")
    log:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/logs_extractTopconsensusedges/cf_{{thr}}.log"
    params:
        analysis=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/",
        cell_order=lambda w: f"{config['datadir']}/{celltype_order_dataset_map[w.dataset]}/inputdata/celltype_order.txt",
        datadir=lambda w: (
            f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
        ),
    shell:
        """
        module load bioconda-2024.02
        Rscript --vanilla {config[scriptdir]}/downstream_network_inference/scMTNI/extractTopconsensusedges_scMTNI.R \
            {params.analysis} cf {wildcards.thr} {params.cell_order} {params.datadir} > {log} 2>&1
        """


rule stabilityKmeansClustering:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.top_{{top}}.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/kmeansclustering_filteredlowexpression/.kmeans_top_{{top}}.done")
    log:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/kmean_log/top_{{top}}.log"
    params:
        path=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis"
    threads: 20
    shell:
        """
        module load matlab-r2023a
        cd {config[scriptdir]}/downstream_network_inference/scMTNI
        matlab -nodisplay -nodesktop -r \
            "path='{params.path}'; prefix='{wildcards.dataset}'; n={wildcards.top}; \
             run stabilityKmeansClustering_scMTNI.m; quit()" \
            > {log} 2>&1
        """

rule gen_bubbleplot_scMTNI:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/kmeansclustering_filteredlowexpression/.kmeans_top_{{top}}.done"
    output:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/kmeansclustering_filteredlowexpression/top{{top}}_bubbles/kmeans_k{{k}}/kmeans_meanprofile_heatmap.pdf"
    params:
        cluster_dir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/kmeansclustering_filteredlowexpression",
        out_dir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/kmeansclustering_filteredlowexpression/top{w.top}_bubbles",
        prefix=lambda w: f"{w.dataset}_stability_kmeans_top{w.top}_k{w.k}"
    conda: scMTNI_plots_conda
    shell:
        """
        mkdir -p {params.out_dir}

        Rscript {config[scriptdir]}/downstream_network_inference/scMTNI/gen_bubbleplot_scMTNI_top.R \
            {wildcards.k} \
            {params.cluster_dir} \
            {params.out_dir} \
            {params.prefix}
        """


###
## LDA
###

rule prepareNetmatrix_scMTNI:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.netfilter.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/.prepare.done")
    params:
        cell_order=lambda w: f"{config['datadir']}/{celltype_order_dataset_map[w.dataset]}/inputdata/celltype_order.txt",
        analysis=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/"
    threads: 5
    shell:
        """
        module load bioconda-2024.02
        Rscript {config[scriptdir]}/downstream_network_inference/scMTNI/prepareNetmatrix_scMTNI_cf_remove_TFA_target.R \
            {params.cell_order} {params.analysis} {wildcards.filter}
        """


rule LDA_analysis_scMTNI:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/.prepare.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/.lda{{filter}}.done")
    log:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/lda{{filter}}.log"
    params:
        cell_order=lambda w: f"{config['datadir']}/{celltype_order_dataset_map[w.dataset]}/inputdata/celltype_order.txt",
        indir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/"
    threads: 20
    shell:
        """
        module load matlab-r2023a
        cd {config[scriptdir]}/downstream_network_inference/scMTNI
        matlab -nodisplay -nodesktop -r \
        "cells=importdata('{params.cell_order}'); indir='{params.indir}'; cf='{wildcards.filter}'; k_lst=[10]; run LDA_analysis_wrapper_scMTNI.m; quit" \
        > {log} 2>&1
        """



rule plotNetworks_LDA:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/.lda{{filter}}.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/.network_all.done")
    params:
        analysis=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis",
        celltype_order=lambda wildcards: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[wildcards.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
        outdir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}/"
    conda: scMTNI_plots_conda
    shell:
        """
        cd {config[scriptdir]}/downstream_network_inference/scMTNI
        Rscript plotNetworks_LDA_v2.R genNetall_multi.R \
            {params.analysis} {params.celltype_order} {wildcards.filter} \
            {wildcards.k} {params.outdir}
        """


rule plotNetworks_LDA_edgetype:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/.lda{{filter}}.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/{{edgetype}}/.network.done")
    params:
        analysis=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis",
        celltype_order=lambda wildcards: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[wildcards.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
        outdir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}/{w.edgetype}",
        datadir=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
    conda: scMTNI_plots_conda
    shell:
        """
        cd {config[scriptdir]}/downstream_network_inference/scMTNI
        Rscript plotNetworks_LDA_v2.R genNetall_multi.R \
            {params.analysis} {params.celltype_order} {wildcards.filter} \
            {wildcards.k} {params.outdir} {params.datadir} {wildcards.edgetype}
        """


rule makeLDATopicRegBubble:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/.network_all.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/.bubble.done")
    params:
        kdir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}",
        celltype_order=lambda wildcards: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[wildcards.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
        datadir=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
    wildcard_constraints: k=r"\d+"
    conda: scMTNI_plots_conda
    shell:
        """
        cd {config[scriptdir]}/downstream_network_inference/scMTNI
        Rscript makeLDATopicRegBubble_ggplot_scMTNI_full_graph.R \
            {params.kdir} {wildcards.k} {params.celltype_order} \
            {params.kdir} {params.datadir}
        """


rule makeLDATopicRegBubble_gene:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/gene/.network.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/gene/.bubble.done")
    params:
        edir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}/gene",
        celltype_order=lambda w: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[w.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
        datadir=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}"
    conda:
        scMTNI_plots_conda
    wildcard_constraints:
        k=r"\d+"
    shell:
        """
        cd {config[scriptdir]}/downstream_network_inference/scMTNI
        Rscript makeLDATopicRegBubble_ggplot_scMTNI_full_graph.R \
            {params.edir} {wildcards.k} {params.celltype_order} \
            {params.edir} {params.datadir}
        """


rule makeLDATopicRegBubble_TFA:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/TFA/.network.done",
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/motifAnalysis/motif_pairwise_pearson_histograms.pdf"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/TFA/.bubble.done")
    params:
        edir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}/TFA",
        celltype_order=lambda w: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[w.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
        datadir=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}",
        motifIndir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/motifAnalysis/",
    conda:
        scMTNI_plots_conda
    wildcard_constraints:
        k=r"\d+"
    shell:
        """
        cd {config[scriptdir]}/downstream_network_inference/scMTNI
        Rscript makeLDATopicRegBubble_ggplot_scMTNI_full_graph.R \
            {params.edir} {wildcards.k} {params.celltype_order} \
            {params.edir} {params.datadir} {params.motifIndir}
        """

# # actually only need to do this for one algo (as expression does not change).
# To identify transcription factors with degenerate motif-derived target profiles, 
# we calculated pairwise Pearson correlations between the weighted TF-to-target profiles 
# within each cell cluster. TF pairs with highly similar profiles (Pearson correlation ≥ 0.995) 
# were represented as edges in an undirected graph, and connected components were defined as motif-degenerate 
# TF families. Because transcription factors with similar DNA-binding motifs may nevertheless 
# differ in their relevance to a given cell state, we further filtered each family using cluster-specific gene expression. 
# TFs classified as lowly expressed in a cluster were removed, and families containing fewer than two remaining TFs were excluded. 
# The resulting families therefore represent groups of expressed TFs whose motif-derived target profiles are effectively
# indistinguishable within the corresponding cell cluster.
rule identifyDegenerateMotif:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/.netfilter.done"
    output:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/motifAnalysis/motif_pairwise_pearson_histograms.pdf"
    conda: scMTNI_plots_conda
    params:
        celltype_order=lambda wildcards: (
            f"{config['datadir']}/"
            f"{celltype_order_dataset_map[wildcards.dataset]}/"
            f"inputdata/celltype_order.txt"
        ),
        indir=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/",
        datadir=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}/",
        suffix = "_network.txt",
        cutoff = 0.995, # may be different if the prior is constructed in a way other than motif
        outdir = lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/motifAnalysis/",
    shell:
        """
        Rscript {config[scriptdir]}/downstream_network_inference/scMTNI/identify_degenerate_motifs.R \
        {params.indir} {params.celltype_order} {params.suffix} {params.cutoff} {params.datadir} {params.outdir}
        """

# enrichment. We ignore the regulators and focus only on the target gene sets.
rule LDA_GO_enrichment_target_only:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/.network_all.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}_target_only/.go_enrich.done")
    wildcard_constraints:
        k=r"\d+"
    params:
        outdir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}_target_only",
        indir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}",
        goterms=lambda w: GoTerms[w.dataset],
        union=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}/Union.txt",
        cell_order=lambda w: f"{config['datadir']}/{celltype_order_dataset_map[w.dataset]}/inputdata/celltype_order.txt",
        enrichAnalyzer = enrichAnalyzer
    shell:
        """
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/lda_topic_net_enrichAnalyzer_wrapper.sh \
            {params.outdir} \
            {wildcards.k} \
            {params.goterms} \
            {params.union} \
            {params.cell_order} \
            T \
            {params.indir} \
            F \
            {params.enrichAnalyzer}
        """

rule LDA_GO_enrichment_edgetype_target_only:
    input:
        f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/{{edgetype}}/.network.done"
    output:
        touch(f"{config['outdir']}/algos/scMTNI/{{algo}}/{{dataset}}/4_1/analysis/lda_TFcellbygene/network{{filter}}/k{{k}}/{{edgetype}}/.go_enrich.done")
    wildcard_constraints:
        k=r"\d+",
        edgetype="gene|TFA"
    params:
        edir=lambda w: f"{config['outdir']}/algos/scMTNI/{w.algo}/{w.dataset}/4_1/analysis/lda_TFcellbygene/network{w.filter}/k{w.k}/{w.edgetype}",
        goterms=lambda w: GoTerms[w.dataset],
        union=lambda w: f"{config['datadir']}/{w.dataset}/inputdata/scMTNI/data/{w.algo}/Union.txt",
        cell_order=lambda w: f"{config['datadir']}/{celltype_order_dataset_map[w.dataset]}/inputdata/celltype_order.txt",
        enrichAnalyzer = enrichAnalyzer
    shell:
        """
        bash {config[scriptdir]}/downstream_network_inference/scMTNI/lda_topic_net_enrichAnalyzer_wrapper.sh \
            {params.edir} \
            {wildcards.k} \
            {params.goterms} \
            {params.union} \
            {params.cell_order} \
            T \
            {params.edir} \
            F \
            {params.enrichAnalyzer}
        """