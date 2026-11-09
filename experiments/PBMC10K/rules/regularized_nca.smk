
rule run_regularized_nca:
    input:
        r=expand("{datadir}/{{dataset}}/all_regulators.txt", datadir=config["datadir"]),
        d=expand("{datadir}/{{dataset}}/normalized_{{c}}_exp.txt", datadir=config["datadir"]),
        g=expand("{datadir}/{{dataset}}/all_genes.txt", datadir=config["datadir"]),
        p=expand("{datadir}/{{dataset}}/pbmc10k_prior.txt", datadir=config["datadir"]),
    output:
        out_A=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}/{{stability_rep}}/adj.txt", outdir=config["outdir"]),
        out_P=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}/{{stability_rep}}/tfa.txt", outdir=config["outdir"]),
        outdir=directory(expand("{outdir}/regularized_nca/{{dataset}}/{{c}}/{{stability_rep}}/", outdir=config["outdir"]))
    params:
        exe="/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/EstimateNCA/NCALearner"
    benchmark: 
        "benchmarks/regularized_nca/{dataset}/{c}/{stability_rep}.benchmark.txt"
    shell:
        """
        {params.exe} -d {input.d} -r {input.r} -g {input.g} -p {input.p} -o {output.outdir} \
        -l {config[algo][regularized_nca][l]}
        """

rule regularized_nca_get_average:
    input:
        in_P=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}/{stability_rep}/tfa.txt", \
            outdir=config["outdir"], stability_rep=["rep_"+str(_) for _ in range(num_random_init)])
    output:
        file_list=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}_stability/tfa_filelist.txt", outdir=config["outdir"]),
        avg_tfa=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}_stability/P.txt", outdir=config["outdir"])
    conda:
        "/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/.conda/envs/python2.7"
    shell:
        """
        echo '{input.in_P}' | tr ' ' '\\n' > {output.file_list}
        python {config[scriptdir]}/regularized_nca/takeAvg.py --tfalist {output.file_list} \
            --tfaout {output.avg_tfa}
        """

rule regularized_nca_get_consensus_net: # this return edge list
    input:
        in_A=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}/{stability_rep}/adj.txt", \
            outdir=config["outdir"], stability_rep=["rep_"+str(_) for _ in range(num_random_init)])
    output:
        file_list=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}_stability/adj_filelist.txt", outdir=config["outdir"]),
        out_adj=expand("{outdir}/regularized_nca/{{dataset}}/{{c}}_stability/A.txt", outdir=config["outdir"])
    params:
        exe="/mnt/dv/wid/projects2/Roy-common/programs/programs/estimateedgeconf/estimateEdgeConf"
    shell:
        """
        echo '{input.in_A}' | tr ' ' '\\n' > {output.file_list}
        {params.exe} {output.file_list} 0.8 {output.out_adj}.tmp alledges
        mv {output.out_adj}.tmpalledge.txt {output.out_adj}
        rm {output.out_adj}.tmp*
        """
