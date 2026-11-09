library(BITFAM)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

in_exp=args[[1]]#"inputtoTFA/0.2/k20/Calb/rep0/E.txt"
tf_lst=args[[2]]#"/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt"
prior_path=args[[3]]#"inputtoTFA/0.2/k20/Calb_prior.txt"
out_A_path=args[[4]]
out_P_path=args[[5]]

# inputtoTFA/0.2/k20/Calb_prior.txt
# testoutput/0.2/Calb/genes.txt 

BITFAM_v2 <- function(data, TF_used, prior_tab, ncores, iter = 8000, tol_rel_obj = 0.005){
    gene_list <- list()                                                                                                                       
    for (i in TF_used) { 
        tmp_gene <- (prior_tab[prior_tab$V1==i,])$V2
        gene_list[[which(TF_used == i)]] <- rownames(data)[rownames(data) %in%                                                                
            tmp_gene]
    }

    X <- t(as.matrix(data))
    chipseq_weight <- matrix(1, nrow = ncol(X), ncol = length(TF_used)) 
    for (i in 1:length(TF_used)) {
        chipseq_weight[, i] <- ifelse(colnames(X) %in% gene_list[[i]],
            1, 0)
    }
    Mask_matrix <- chipseq_weight

    N <- dim(X)[1]
    D <- dim(X)[2]
    K <- length(TF_used)
    data_to_model <- list(N = N, D = D, K = K, X = X, Mask = Mask_matrix) 
    library(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = ncores)
    set.seed(100)
    pca_beta_piror <- "\ndata {\nint<lower=0> N; // Number of samples\nint<lower=0> D; // The original dimension\nint<lower=0> K; // The latent dimension\nmatrix[N, D] X; // The data matrix\nmatrix[D, K] Mask; // The binary mask of prior knowledge indicate the target of TFs\n}\n\nparameters {\nmatrix<lower=0, upper=1>[N, K] Z; // The latent matrix\nmatrix[D, K] W; // The weight matrix\nreal<lower=0> tau; // Noise term\nvector<lower=0>[K] alpha; // ARD prior\n}\n\ntransformed parameters{\nmatrix<lower=0>[D, K] t_alpha;\nreal<lower=0> t_tau;\nfor(wmd in 1:D){\nfor(wmk in 1:K){\nt_alpha[wmd, wmk] = Mask[wmd, wmk] == 1 ? inv(sqrt(alpha[wmk])) : 0.01;\n}\n}\nt_tau = inv(sqrt(tau));\n}\nmodel {\ntau ~ gamma(1,1);\nto_vector(Z) ~ beta(0.5, 0.5);\nalpha ~ gamma(1e-3,1e-3);\nfor(d in 1:D){\nfor(k in 1:K){\nW[d,k] ~ normal(0, t_alpha[d, k]);\n}\n}\nto_vector(X) ~ normal(to_vector(Z*W'), t_tau);\n} "
    m_beta_prior <- stan_model(model_code = pca_beta_piror)
    suppressWarnings(fit.vb <- vb(m_beta_prior, data = data_to_model,
        algorithm = "meanfield", iter = iter, output_samples = 300,
        tol_rel_obj = tol_rel_obj))
    BITFAM_list <- list(Model = fit.vb, TF_used = TF_used, Genes = rownames(data),                                                           
        Cell_names = colnames(data))
    return(BITFAM_list)
}

# require feature by cell
raw_dat = read.table(in_exp, check.names=FALSE, sep='\t', row.names=1)
TF_used = read.table(tf_lst, check.names=FALSE)[,"V1"]
prior_tab = read.table(prior_path, check.names=FALSE, sep='\t')
colnames(raw_dat) <- paste0("cell", 1:ncol(raw_dat))

# here we just replicate BITFAM_preprocess as it errs and we don't want filtering
process_data <- CreateSeuratObject(counts = raw_dat) # error occurs in this step, so I omitted min.cells and min.features
process_data <- NormalizeData(object = process_data) 
process_data <- FindVariableFeatures(object = process_data, nfeatures = 5000) # we have way fewer features so this step does nothing
data_matrix_normalized <- as.matrix(GetAssayData(object = process_data)[VariableFeatures(process_data), 
        ])
# note that the resulting data_matrix_normalized has rows (genes) reordered. Need to address that
rownames(data_matrix_normalized) <- VariableFeatures(process_data)
colnames(data_matrix_normalized) <- colnames(GetAssayData(object = process_data))

prior_tab <- prior_tab[prior_tab$V2 %in% rownames(data_matrix_normalized),]
# Below, I just follow the BITFAM function
BITFAM_res <- BITFAM_v2(data_matrix_normalized, TF_used, prior_tab, parallel::detectCores(), iter = 8000, tol_rel_obj = 0.005)

Z <- BITFAM_activities(BITFAM_res)

W <- apply(extract(BITFAM_res$Model, "W")[[1]], c(2,3), mean)
rownames(W) <- BITFAM_res$Genes
colnames(W) <- BITFAM_res$TF_used

Z_ordered <- t(Z)
Z_ordered <- Z_ordered[TF_used, colnames(raw_dat)]

W_ordered <- W[rownames(raw_dat), TF_used]

write.table(format(W_ordered, digits=5), out_A_path, quote = FALSE, sep = "\t",row.names = FALSE,
            col.names = FALSE)
write.table(format(Z_ordered, digits=5), out_P_path, quote = FALSE, sep = "\t",row.names = FALSE,
            col.names = FALSE)