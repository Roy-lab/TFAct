# use k means to control the number of clusters

library(Seurat)
library(dplyr)
library(ggplot2)
library(mclust)

args = commandArgs(trailingOnly=TRUE)
inferred_P = args[[1]]#TFA_output/baselines/hals_prior/0.2/k20/Anc2/rep0/E_stability/P.txt
ground_truth_cell_type = args[[2]]#inputtoTFA/0.2/k20/Anc2/rep0/cell_celltype_assignments.txt
tf_lst = args[[3]] #/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt
P_sep = args[[4]] # seperator for inferred TFAs. can pass ' ' or  $'\t'
rownames_first_col = args[[5]] # if T: read the first column as row names. Else there are no row names.
out_stat = args[[6]]
# exp_mat = args[[5]] #inputtoTFA/0.2/k20/Anc2/rep0/E.txt
thresholding_noises = "F"
if(length(args) > 6){
    thresholding_noises <- args[[7]] 
}

npca_c = 10 # as we have a few features. This is for TFA's pca, but note that we run kmeans using all TFA features


ct = read.table(ground_truth_cell_type, check.names=FALSE, header=FALSE)
tf = read.table(tf_lst, check.names=FALSE, header=FALSE)[[1]]
if(rownames_first_col=="T"){
    P = read.table(inferred_P, row.names=1, check.names=FALSE, header=FALSE, sep = P_sep)
    P = P[tf, ]
} else {
    P = read.table(inferred_P, check.names=FALSE, header=FALSE, sep = P_sep)
}

if(thresholding_noises=="T"){
    P = P / rowSums(abs(P))
    P[P<1/ncol(P)]=0
}


all_celltypes = unique(ct)

rownames(ct) <- paste0("cell",c(0:(ncol(P)-1)))
colnames(ct) <- c("cellType")

colnames(P) <- paste0("cell",c(0:(ncol(P)-1)))
rownames(P) <- tf


P <- CreateSeuratObject(P, assay="TFA", meta.data = ct, project = "simulation")

# P <- NormalizeData(P, normalization.method = "LogNormalize", scale.factor = 10000)

P <- FindVariableFeatures(P, selection.method = "vst", nfeatures = 2000) # this does nothing as we don't have that many features

# scale data
all.genes <- rownames(P)
P <- ScaleData(P, features = all.genes)

P <- RunPCA(P, features = VariableFeatures(object = P), npcs=npca_c, nfeatures.print=npca_c)

# P <- FindNeighbors(P, features = rownames(P@assays$TFA@data)) # we can't specify # clusts
# P <- FindClusters(P, resolution = 0.5)

set.seed(123)
km.res <- kmeans(t(as.matrix(P@assays$TFA@data)), 10, nstart = 25)

ARI <- adjustedRandIndex(km.res[["cluster"]], ct[[1]])

outdf <- data.frame(ARI=c(ARI))

write.table(outdf, out_stat, quote = FALSE, sep = "\t",row.names = FALSE,
            col.names = TRUE)

# TODO
# P@meta.data$inferredCellType=km.res[["cluster"]]


# # for visualization
# P <- RunUMAP(P, dims = 1:10, min.dist = 0.8, spread = 1.5, n.neighbors=50)

# # ground truth
# DimPlot(P, reduction = "umap", group.by="cellType")
# ggsave("test2_pca.png")

# DimPlot(P, reduction = "umap", group.by="inferredCellType")
# ggsave("test3_pca.png")





# # now using the original expression data set
# exp = read.table(exp_mat, row.names=1, check.names=FALSE, header=FALSE)
# colnames(exp) <- rownames(ct)
# exp_npca_c = length(tf) # As we perform kmeans on pca, we set the number of features to the number of tfs, consistent with TFA dimension



# E <- CreateSeuratObject(exp, assay="RNA", meta.data = ct, project = "simulatedExpression")
# E <- NormalizeData(E, normalization.method = "LogNormalize", scale.factor = 10000)
# E <- FindVariableFeatures(E, selection.method = "vst", nfeatures = 2000) # this does nothing as we don't have that many features
# E <- ScaleData(E, features = rownames(E))

# E <- RunPCA(E, features = VariableFeatures(object = E), npcs=exp_npca_c, nfeatures.print=10)

# km.res.exp <- kmeans(E@reductions$pca@cell.embeddings, 10, nstart = 25)

# ARI_exp <- adjustedRandIndex(km.res.exp[["cluster"]], ct[[1]])