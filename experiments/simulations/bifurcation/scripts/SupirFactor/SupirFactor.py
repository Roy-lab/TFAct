import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_exp", type=str, required=True)
    parser.add_argument("--tf_lst", type=str, required=True)
    parser.add_argument("--prior_path", type=str, required=True)
    parser.add_argument("--out_A", type=str, required=True)
    parser.add_argument("--out_P", type=str, required=True)
    parser.add_argument("--model_used", help="Hierarchical or Shallow",type=str, default="Hierarchical")

    args = parser.parse_args()



    from supirfactor.setup import dataloader, prior as setup_prior, run as setup_run
    import supirfactor.run as sf_run
    from supirfactor import erv
    import supirfactor.postprocess as sf_post # not working for some scanpy version
    import pandas as pd
    import anndata as ad
    import numpy as np
    import supirfactor.models as md


    data_file=args.in_exp#"/mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/simulation/072024_cell_lineage_formulation_w_lda_noise/inputtoTFA/0.2/k20/Calb/rep0/E.txt"
    tf_lst=args.tf_lst#"/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt"
    prior_path=args.prior_path
    model_used=args.model_used

    adata = pd.read_csv(data_file, sep='\t', index_col=0, header=None).astype(np.float32)
    adata.columns = ["cell"+str(i) for i in range(adata.shape[1])]
    adata = adata.transpose()# need cell by gene
    adata = ad.AnnData(adata, dtype=adata.values.dtype)

    prior_tab = np.loadtxt(prior_path, dtype=str)
    genes=adata.var_names
    tfs=np.loadtxt(tf_lst, dtype=str)

    prior = np.zeros((len(genes), len(tfs)))

    gene_dict={genes[i]:i for i in range(len(genes))}
    reg_dict={tfs[i]:i for i in range(len(tfs))}

    for tf, g, _ in prior_tab:
        prior[gene_dict[g], reg_dict[tf]] = _

    prior= pd.DataFrame(prior, index=genes, columns=tfs)
    prior = prior.loc[:, (prior.sum(0) != 0)]

    omitted_tfs = np.setdiff1d(tfs, prior.columns)

    cv_prior = False

    use_prior, use_gs, use_prior_map = setup_prior.split(prior, cv_prior=cv_prior, fraction_gs=0.2)

    # well this step is not working
    if model_used=="Hierarchical":
        model = md.Hierarchical(use_prior, pi_structure=[len(tfs)]) # from the paper
    elif model_used=="Shallow":
        model = md.Shallow(use_prior, pi_structure=[len(tfs)]) 
    else:
        print("--model_used should be either Hierarchical or Shallow")
        exit()

    run_configs = setup_run.defaults(model)
    sample_sets = dataloader.split_samples(adata, validation_size=0.1) # not validating. All training

    randomize_data = False
    preprocessor=None
    loaders = dataloader.setup_loaders(adata, batch_size=512, **sample_sets, randomize_data=randomize_data, verbose=True)

    tr_loss = sf_run.fit(model, loaders[0], epochs=200, validation_loader=loaders[1], verbose=True, **run_configs)


    # tr_loss = sf_post.loss2df(tr_loss)
    # __ = tr_loss.plot.scatter(x='epoch', y='R2')
    # __ = tr_loss.plot.scatter(x='epoch', y='trl')
    # __ = tr_loss.plot.scatter(x='epoch', y='vl')


    model_erv = erv.ERV(model, loaders[1]) # representation of GRN

    ##
    # Get GRN
    ##
    grn = model_erv['0_erv_df']

    grn.columns = use_prior.columns
    grn.index = adata.var_names


    ##
    # Compute TFA
    ##

    LA = sf_post.generate_latent_embedding(model, adata, batch_size=512, return_nodes={'W': 'TFA'})
    TFA = sf_post.to_anndata(LA['TFA'], var_names=use_prior.columns.tolist(), obsm=adata.obsm, uns=adata.uns, obs=adata.obs)
    TFA = pd.DataFrame(TFA.X,index=TFA.obs_names, columns=TFA.var_names)
    TFA = TFA.transpose()

    out_TFA = pd.DataFrame(np.zeros((len(tfs), adata.X.shape[0])), index=tfs, columns=adata.obs_names)
    out_TFA.loc[TFA.index, TFA.columns] = TFA # in case there are missing tfs, we set their tfa to 0

    out_grn = pd.DataFrame(np.zeros((adata.X.shape[1], len(tfs))), index=genes, columns=tfs)
    out_grn.loc[grn.index, grn.columns] = grn

    np.savetxt(args.out_A, out_grn.to_numpy(), fmt='%1.5f', delimiter='\t')
    np.savetxt(args.out_P, out_TFA.to_numpy(), fmt='%1.5f', delimiter='\t')
