# full slingshot for pseudotime estimation and tradeseq analysis. Same seed as TIGER_zenodo/code/PBMC10K_v2.R
library(Seurat)
library(dplyr)
library(ggplot2)
library(OmnipathR)
library(slingshot)
library(parallel)
# BiocManager::install("tradeSeq") # run this if not done
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(presto)
library(igraph)


max_cores = 50

set.seed(1234) 

outdir="TFA_output"

celltypes=c("CD4_Naive","CD4_TCM","CD8_Naive","CD16_Mono","NK","Treg","CD14_Mono","cDC","CD8_TEM_1","Intermediate_B",
        "Memory_B","Plasma","CD4_TEM","MAIT","Naive_B","CD8_TEM_2","pDC","gdT","HSPC")
cellnames = lapply(celltypes, function(x){read.table(paste0('TIGER_zenodo/data/PBMC_ct/', x, '_cellnames.txt'))$V1})
names(cellnames) <- celltypes
stability_selection_algos=c('HALS-p', 'ibpg_2fct', 'ibpg', 'ibpg_2fct_500_1000', 'ibpg_2fct_500_0')
algo_w_TF_names=c('regularized_nca')
other_algos=c('TIGER','supirfactor_hierarchical', 'supirfactor_shallow', 'NetREX_CF','NINCA_quad') 
algos = c(stability_selection_algos, algo_w_TF_names, other_algos)

Ps_cellnames <- unlist(cellnames) # TFA cell names follow this order



###
# Step 1: Check that slingshot does not produce biological plausible trajectories using the full data
###
# Create Seurat object. from TIGER's zenodo dir
pbmc <- readRDS("TIGER_zenodo/data/PBMC/pbmc_final.rds")

dimred <- pbmc@reductions$wnn.umap@cell.embeddings
clustering <- pbmc$celltype

sds_umap <- slingshot(dimred, clusterLabels = clustering, start.clus = "HSPC")

pca_emb <- pbmc@reductions$pca@cell.embeddings
clustering <- pbmc$celltype

sds_pca_list <- list()

for (npc in c(3, 5)) {
  sds_pca_list[[paste0("PC", npc)]] <- slingshot(
    pca_emb[, 1:npc, drop=FALSE],
    clusterLabels = clustering,
    start.clus = "HSPC"
  )

  cat("\n", npc, "PCs\n")
  print(slingLineages(sds_pca_list[[paste0("PC", npc)]]))
}


pdf("slingshot_trees.pdf", width=10, height=7, useDingbats=FALSE)
par(family="Helvetica", mar=c(1,1,3,1))

plot_one_tree <- function(sds, main, root="HSPC") {
  lineage_dat <- slingLineages(sds)

  edge_list <- do.call(rbind, lapply(lineage_dat, function(l) {
    cbind(parent=l[-length(l)], child=l[-1])
  }))

  edge_df <- unique(as.data.frame(edge_list, stringsAsFactors=FALSE))
  g <- graph_from_data_frame(edge_df, directed=TRUE)

  lay <- layout_as_tree(g, root=which(V(g)$name == root))
  lay[,1] <- lay[,1] * 2.5
  lay[,2] <- lay[,2] * 1.5

  # Reverse vertical direction so root is on top
  lay[,2] <- -lay[,2]

  labels <- V(g)$name
  box_w <- 0.18 * nchar(labels) + 0.7
  box_h <- 0.55

  xlim <- range(lay[,1] + cbind(-box_w/2, box_w/2)) + c(-0.5, 0.5)
  ylim <- range(lay[,2]) + c(-1, 1)

  plot(NA, xlim=xlim, ylim=ylim, axes=FALSE, xlab="", ylab="", main=main)

  # edges first
  for (i in seq_len(nrow(edge_df))) {
    p <- match(edge_df$parent[i], labels)
    c <- match(edge_df$child[i], labels)

    x1 <- lay[p,1]
    y1 <- lay[p,2] - box_h/2
    x2 <- lay[c,1]
    y2 <- lay[c,2] + box_h/2

    arrows(x1, y1, x2, y2, length=0.08, lwd=1.5, col="grey50")
  }

  # boxes + labels
  for (i in seq_along(labels)) {
    rect(
      lay[i,1] - box_w[i]/2,
      lay[i,2] - box_h/2,
      lay[i,1] + box_w[i]/2,
      lay[i,2] + box_h/2,
      col="lightblue",
      border="black",
      lwd=1.2
    )

    text(lay[i,1], lay[i,2], labels[i], cex=0.9)
  }
}

for (nm in names(sds_pca_list)) {
  plot_one_tree(sds_pca_list[[nm]], paste("Slingshot tree:", nm))
}

plot_one_tree(sds_umap, "Slingshot tree: WNN UMAP")

dev.off()


##
# Step 2: As slingshot on the original data give biologically implausible lineages, we define
# lineages separately and then apply slingshot
##
lineage_sets <- list(
  Monocyte = list(
    cells = c("HSPC", "CD14 Mono", "CD16 Mono"),
    end = "CD16 Mono"
  ),
  CD4 = list(
    cells = c("HSPC", "CD4 Naive", "CD4 TCM", "CD4 TEM"),
    end = "CD4 TEM"
  ),
  CD8 = list(
    cells = c("HSPC", "CD8 Naive", "CD8 TEM_1", "CD8 TEM_2"),
    end = NULL
  ),
  B = list(
    cells = c("HSPC", "Naive B", "Intermediate B", "Memory B"),
    end = "Memory B"
  )
)

sds_list <- list()
pt_list <- list()
cw_list <- list()

# inferred lineages should match lineage_sets 
for (nm in names(lineage_sets)) {

  x <- lineage_sets[[nm]]
  keep <- clustering %in% x$cells

  if (is.null(x$end)) {
    sds <- slingshot(
      dimred[keep, , drop=FALSE],
      clusterLabels = clustering[keep],
      start.clus = "HSPC"
    )
  } else {
    sds <- slingshot(
      dimred[keep, , drop=FALSE],
      clusterLabels = clustering[keep],
      start.clus = "HSPC",
      end.clus = x$end
    )
  }

  sds_list[[nm]] <- sds
  pt_list[[nm]] <- slingPseudotime(sds)
  cw_list[[nm]] <- slingCurveWeights(sds)

  cat("\n", nm, "\n")
  print(slingLineages(sds))
}

# check pseudotime
pdf("slingshot_pseudotime_by_celltype.pdf", width=8, height=5)

for (nm in names(pt_list)) {
  pt <- pt_list[[nm]][,1]
  cells <- rownames(pt_list[[nm]])
  ct <- clustering[cells]

  boxplot(
    pt ~ ct,
    las=2,
    main=nm,
    ylab="Slingshot pseudotime",
    xlab=""
  )
}
dev.off()

plot_slingshot_tree <- function(sds, main, root="HSPC") {

  lineages <- slingLineages(sds)

  edges <- unique(do.call(rbind, lapply(lineages, function(l) {
    data.frame(from=head(l, -1), to=tail(l, -1))
  })))

  g <- graph_from_data_frame(edges, directed=TRUE)

  lay <- layout_as_tree(g, root=which(V(g)$name == root))
  lay[, 2] <- -lay[, 2]

  plot(
    g,
    layout=lay,
    main=main,
    vertex.shape="rectangle",   # Rectangular nodes
    vertex.size=95,             # Rectangle width
    vertex.size2=18,            # Rectangle height
    vertex.color="lightblue",
    vertex.frame.color="black",
    vertex.label.color="black",
    vertex.label.cex=1.5,
    vertex.label.dist=0,
    edge.arrow.size=2,
    edge.color="grey50",
    edge.width=4,
    margin=0.1
    )
}

# Plot all four lineage trees in one PDF
dir.create("plots/tradeseq", recursive=TRUE, showWarnings=FALSE)

pdf(
  "plots/tradeseq/slingshot_lineage_trees.pdf",
  width=14, height=10,
  family="Helvetica"
)

par(mfrow=c(2, 2), mar=c(2, 2, 3, 2))

for (nm in names(sds_list)) {
  plot_slingshot_tree(sds_list[[nm]], main=nm)
}

dev.off()
#
## Last check: see if the slingshot works on the combined clusters by specifying end clusters
#
all_lineage_cells <- unique(unlist(lapply(lineage_sets, function(x) x$cells)))

keep <- clustering %in% all_lineage_cells

end_clusters <- c(
  "CD16 Mono",
  "CD4 TEM",
  "CD8 TEM_2",
  "Memory B"
)

sds_combined <- slingshot(
  dimred[keep, , drop=FALSE],
  clusterLabels = clustering[keep],
  start.clus = "HSPC",
  end.clus = end_clusters
)

slingLineages(sds_combined) # it still produced inconsistent lineages 
# pca
for (npc in c(3, 5)) {
  pca_red <- pbmc@reductions$pca@cell.embeddings[keep, 1:npc, drop=FALSE]
  sds_pca <- slingshot(
    pca_red,
    clusterLabels=clustering[keep],
    start.clus="HSPC",
    end.clus=end_clusters
  )
  cat("\n###", npc, "PCs ###\n")
  print(slingLineages(sds_pca)) # still wrong
}

###
# Step 3. Tradeseq analysis. We used the results of applying slingshot separately to each lineage.
###

RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))

all_reg <- read.table("TIGER_zenodo/data/PBMC_ct/all_regulators.txt")$V1

# per vignette requires counts
counts_all <- GetAssayData(pbmc, assay="RNA", layer="counts")


##
## Tradeseq 1. Build the joint pseudotime and weight matrices
##
library(mgcv)
library(Matrix)

lineages <- names(pt_list)

all_cells <- unique(unlist(lapply(pt_list, rownames)))

# preserve PBMC cell order. 9171 cells out of the 10412 cells 
all_cells <- colnames(pbmc)[colnames(pbmc) %in% all_cells]

pt_joint <- matrix(
  0,
  nrow=length(all_cells),
  ncol=length(lineages),
  dimnames=list(all_cells, lineages)
)

cw_joint <- matrix(
  0,
  nrow=length(all_cells),
  ncol=length(lineages),
  dimnames=list(all_cells, lineages)
)

for (nm in lineages) {

  cells <- rownames(pt_list[[nm]])

  pt <- pt_list[[nm]][,1]
  cw <- cw_list[[nm]][,1]

  # normalize pseudotime to [0,1]
  pt <- (pt - min(pt, na.rm=TRUE)) /
        (max(pt, na.rm=TRUE) - min(pt, na.rm=TRUE))

  pt_joint[cells,nm] <- pt
  cw_joint[cells,nm] <- cw
}

# Important: HSPCs occur in several independently fitted trajectories.
# Normalize so total lineage membership for each cell is 1.
cw_joint <- cw_joint / rowSums(cw_joint)

###
## Tradeseq 2. TF expression, within lineage
###
counts <- counts_all[, all_cells, drop=FALSE]

tf_genes <- intersect(all_reg, rownames(counts)) # 335
tf_idx <- match(tf_genes, rownames(counts))

set.seed(1234)

# ～ 2 hrs
icMat <- evaluateK(
  counts=counts,
  pseudotime=pt_joint,
  cellWeights=cw_joint,
  k=3:7,
  nGenes=200,
  verbose=TRUE
)

mean_aic <- colMeans(icMat, na.rm=TRUE) 
mean_aic # best: k = 6
#     k: 3     k: 4     k: 5     k: 6     k: 7 
# 3385.164 3383.640 3376.877 3374.340 3389.450 

set.seed(1234)
numknots <- 6

rna_fit <- fitGAM(
  counts=counts,
  genes=tf_idx,
  pseudotime=pt_joint,
  cellWeights=cw_joint,
  nknots=numknots,
  verbose=TRUE
)

table(rowData(rna_fit)$tradeSeq$converged) # 321 TFs converged
# FALSE  TRUE 
#    14   321 

# within lineage

# See the doc
rna_assoc <- associationTest(rna_fit, lineage=TRUE)
rna_assoc <- as.data.frame(rna_assoc)

rna_start_end <- startVsEndTest(
  rna_fit,
  global=TRUE,
  lineages=TRUE
)


dir.create("plots/tradeseq/DE/startVsEnd", recursive=TRUE, showWarnings=FALSE)
dir.create("plots/tradeseq/DE/association", recursive=TRUE, showWarnings=FALSE)

write.table(
  rna_start_end,
  "plots/tradeseq/DE/startVsEnd/results.txt",
  sep="\t", quote=FALSE, col.names=NA
)

write.table(
  rna_assoc,
  "plots/tradeseq/DE/association/results.txt",
  sep="\t", quote=FALSE, col.names=NA
)

# overall test, start vs end
oStart <- order(rna_start_end$waldStat, decreasing=TRUE)

for (i in 1:10) {
  gene <- names(rna_fit)[oStart[i]]
  pdf(paste0("plots/tradeseq/DE/startVsEnd/", i, "_", gene, ".pdf"), 7, 5)
  print(plotSmoothers(rna_fit, counts, gene=gene))
  dev.off()
}

# lineage-specific, top 5.
# lineage 1-4 are found as colnames(pt_joint) # [1] "Monocyte" "CD4" "CD8" "B"
for (j in 1:4) {
  dir.create(
    paste0("plots/tradeseq/DE/startVsEnd/lineage", j),
    recursive=TRUE, showWarnings=FALSE
  )
  o <- order(rna_start_end[[paste0("waldStat_lineage", j)]], decreasing=TRUE)
  for (i in 1:5) {
    gene <- names(rna_fit)[o[i]]
    pdf(
      paste0(
        "plots/tradeseq/DE/startVsEnd/lineage", j, "/",
        i, "_", gene, ".pdf"
      ),
      7, 5
    )
    print(plotSmoothers(rna_fit, counts, gene=gene))
    dev.off()
  }
}

## save top 10 startvsend, overall and lineage specific
lineage_names <- colnames(pt_joint)

top10 <- data.frame(
  overall=names(rna_fit)[order(rna_start_end$waldStat, decreasing=TRUE)[1:10]]
)

for (j in seq_along(lineage_names)) {
  o <- order(rna_start_end[[paste0("waldStat_lineage", j)]], decreasing=TRUE)
  top10[[lineage_names[j]]] <- names(rna_fit)[o[1:10]]
}

write.table(
  top10,
  "plots/tradeseq/DE/startVsEnd/top10_genes.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

# Overall, association
oAssoc <- order(rna_assoc$waldStat, decreasing=TRUE)

dir.create("plots/tradeseq/DE/association", recursive=TRUE, showWarnings=FALSE)

for (i in 1:10) {
  gene <- names(rna_fit)[oAssoc[i]]
  pdf(paste0("plots/tradeseq/DE/association/", i, "_", gene, ".pdf"), 7, 5)
  print(plotSmoothers(rna_fit, counts, gene=gene))
  dev.off()
}

# Lineage-specific top 5 plots
lineage_names <- colnames(pt_joint)

for (j in seq_along(lineage_names)) {
  dir.create(
    paste0("plots/tradeseq/DE/association/lineage", j),
    recursive=TRUE, showWarnings=FALSE
  )

  o <- order(rna_assoc[[paste0("waldStat_", j)]], decreasing=TRUE)

  for (i in 1:5) {
    gene <- names(rna_fit)[o[i]]
    pdf(
      paste0(
        "plots/tradeseq/DE/association/lineage", j, "/",
        i, "_", gene, ".pdf"
      ),
      7, 5
    )
    print(plotSmoothers(rna_fit, counts, gene=gene))
    dev.off()
  }
}

# Save top 10 genes: overall and lineage-specific
top10 <- data.frame(
  overall=names(rna_fit)[order(rna_assoc$waldStat, decreasing=TRUE)[1:10]]
)

for (j in seq_along(lineage_names)) {
  o <- order(rna_assoc[[paste0("waldStat_", j)]], decreasing=TRUE)
  top10[[lineage_names[j]]] <- names(rna_fit)[o[1:10]]
}

write.table(
  top10,
  "plots/tradeseq/DE/association/top10_genes.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

###
## Tradeseq 3. TF expression, between lineages
###
rna_pattern <- patternTest(rna_fit, global=TRUE, pairwise=TRUE)
rna_pattern <- as.data.frame(rna_pattern)

rna_end <- diffEndTest(rna_fit, global=TRUE, pairwise=TRUE)
rna_end <- as.data.frame(rna_end)

dir.create("plots/tradeseq/DE/pattern", recursive=TRUE, showWarnings=FALSE)
dir.create("plots/tradeseq/DE/diffEnd", recursive=TRUE, showWarnings=FALSE)

write.table(rna_pattern, "plots/tradeseq/DE/pattern/results.txt", sep="\t", quote=FALSE, col.names=NA)
write.table(rna_end, "plots/tradeseq/DE/diffEnd/results.txt", sep="\t", quote=FALSE, col.names=NA)
## ---------------------------
## patternTest: top 10 + plots
## ---------------------------
oPattern <- order(rna_pattern$waldStat, decreasing=TRUE)
top10_pattern <- names(rna_fit)[oPattern[1:10]]

for (i in 1:5) {
  gene <- top10_pattern[i]
  pdf(paste0("plots/tradeseq/DE/pattern/", i, "_", gene, ".pdf"), 7, 5)
  print(plotSmoothers(rna_fit, counts, gene=gene))
  dev.off()
}

## ---------------------------
## diffEndTest: top 10 + plots
## ---------------------------
oEnd <- order(rna_end$waldStat, decreasing=TRUE)
top10_end <- names(rna_fit)[oEnd[1:10]]

for (i in 1:5) {
  gene <- top10_end[i]
  pdf(paste0("plots/tradeseq/DE/diffEnd/", i, "_", gene, ".pdf"), 7, 5)
  print(plotSmoothers(rna_fit, counts, gene=gene))
  dev.off()
}

rename_lineage_cols <- function(df, lineage_names) {
  for (i in seq_along(lineage_names)) {
    for (j in seq_along(lineage_names)) {
      old <- paste0(i, "vs", j)
      new <- paste0(lineage_names[i], "_vs_", lineage_names[j])
      colnames(df) <- gsub(old, new, colnames(df), fixed=TRUE)
    }
  }
  df
}

lineage_names <- colnames(pt_joint)

rna_pattern_named <- rename_lineage_cols(rna_pattern, lineage_names)
rna_end_named <- rename_lineage_cols(rna_end, lineage_names)

# top 10 genes
pairs <- combn(seq_along(lineage_names), 2)

top10_pattern <- data.frame(
  overall = rownames(rna_pattern_named)[
    order(rna_pattern_named$waldStat, decreasing=TRUE)[1:10]
  ]
)

top10_end <- data.frame(
  overall = rownames(rna_end_named)[
    order(rna_end_named$waldStat, decreasing=TRUE)[1:10]
  ]
)

for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]
  j <- pairs[2, k]
  nm <- paste0(lineage_names[i], "_vs_", lineage_names[j])

  top10_pattern[[nm]] <- rownames(rna_pattern_named)[
    order(rna_pattern_named[[paste0("waldStat_", nm)]], decreasing=TRUE)[1:10]
  ]

  top10_end[[nm]] <- rownames(rna_end_named)[
    order(rna_end_named[[paste0("waldStat_", nm)]], decreasing=TRUE)[1:10]
  ]
}

write.table(
  top10_pattern,
  "plots/tradeseq/DE/pattern/top10_genes.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

write.table(
  top10_end,
  "plots/tradeseq/DE/diffEnd/top10_genes.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

## ----------------------------------
## transient genes: vignette approach
## ----------------------------------
dir.create("plots/tradeseq/DE/transient", recursive=TRUE, showWarnings=FALSE)

patternRes <- rna_pattern
endRes <- rna_end

patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes <- patternRes[, c("Gene", "pattern")]

endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
endRes <- endRes[, c("Gene", "end")]

compare <- merge(patternRes, endRes, by="Gene", all=FALSE)
compare$transientScore <-
  rank(-compare$end, ties.method="min")^2 + rank(compare$pattern, ties.method="random")^2

compare <- compare[order(compare$transientScore, decreasing=TRUE), ]

write.table(
  compare,
  "plots/tradeseq/DE/transient/results.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

write.table(
  compare[1:10, ],
  "plots/tradeseq/DE/transient/top10_genes.txt",
  sep="\t", quote=FALSE, row.names=FALSE
)

p <- ggplot(compare, aes(x=log(pattern), y=log(end))) +
  geom_point(aes(col=transientScore)) +
  labs(x="patternTest Wald Statistic (log scale)",
    y="diffEndTest Wald Statistic (log scale)") +
  scale_color_continuous(low="yellow", high="red") +
  theme_classic()

ggsave(
  "plots/tradeseq/DE/transient/transient_score.pdf",
  p, width=7, height=5
)

for (i in 1:5) {
  gene <- compare$Gene[i]
  pdf(paste0("plots/tradeseq/DE/transient/", i, "_", gene, ".pdf"), 7, 5)
  print(plotSmoothers(rna_fit, counts, gene=gene))
  dev.off()
}

###
## Tradeseq 4. prepare TFA
###
get_concate_P <- function(alg) {
    if (alg %in% stability_selection_algos){
        cur_P <- do.call(cbind, mclapply(celltypes, function(c){read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'_stability/P.txt'), check.names=FALSE, header=FALSE)}, mc.cores=max_cores))
    } else if (alg %in% algo_w_TF_names){
        cur_P <- do.call(cbind, mclapply(celltypes, function(c){
            cur_P_ <- read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'_stability/P.txt'), check.names=FALSE, header=FALSE)
            rownames(cur_P_) <- cur_P_$V1
            cur_P_$V1 <- NULL
            nonactive_tfs <- setdiff(all_reg, rownames(cur_P_))
            nonactive_tfs <- as.data.frame(matrix(0, length(nonactive_tfs), ncol(cur_P_)), row.names=nonactive_tfs)
            colnames(cur_P_) <- colnames(nonactive_tfs) 
            cur_P_ <- rbind(cur_P_, nonactive_tfs)
            cur_P_ <- cur_P_[all_reg,]
            }, mc.cores=max_cores))
    }
    else if(alg %in% other_algos){
        cur_P <- do.call(cbind, mclapply(celltypes, function(c){read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'/P.txt'), check.names=FALSE, header=FALSE)}, mc.cores=max_cores))
    } else {
        throw("unsupport algorithm name")
    }

    # does not scale and remove low TFA here, as TF with high scaled activity may be low in a global sense
    return(cur_P)
}

get_TFA_matrix <- function(alg) {
  P <- as.matrix(get_concate_P(alg))
  if (is.null(rownames(P)) || !all(all_reg %in% rownames(P))) {
    if (nrow(P) == length(all_reg)) {
      rownames(P) <- all_reg
    } else {
      warning("Row number does not match length(all_reg); keeping existing rownames.")
    }
  }

  if (ncol(P) == length(Ps_cellnames)) {
    colnames(P) <- Ps_cellnames
  } else {
    warning("Column number does not match length(Ps_cellnames); keeping existing colnames.")
  }
  # subset only cells that are actually available
  keep <- intersect(all_cells, colnames(P))
  if (length(keep) < length(all_cells)) {
    warning(length(all_cells) - length(keep), " requested cells are missing.")
  }
  P[, keep, drop=FALSE]
}

formal_algo_names <- function(dat, column_name){
    dat[[column_name]][dat[[column_name]] == "ibpg"] <- "IBPG-p"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct"] <- "TFAct-F_500_100"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_1000"] <- "TFAct-F_500_1000"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_0"] <- "TFAct-F_500_0"
    return(dat)
}

## the source code from https://github.com/statOmics/tradeSeq/blob/master/R/fitGAM.R says that
#' Note, that if modeling non-count data (e.g. setting \code{family = "gaussian"}),
#' you may want to avoid the default TMM normalization by overriding it using
#' the \code{offset} argument.

set.seed(1234) 

TFA_fits <- list()
TFA_assoc <- list()
TFA_start_end <- list()
TFA_pattern <- list()
TFA_end <- list()

nonnegative_algos <- algos[sapply(algos, function(alg) {
  P <- get_TFA_matrix(alg)
  all(P >= 0, na.rm=TRUE)
})]

nonnegative_algos # nonnegative TFA values
# as expected, 
# [1] "HALS-p"                   "ibpg_2fct"               
# [3] "ibpg"                     "ibpg_2fct_500_1000"      
# [5] "ibpg_2fct_500_0"          "TIGER"                   
# [7] "supirfactor_hierarchical"

for (alg in nonnegative_algos) {
  cat("\nRunning:", alg, "\n")
  P <- get_TFA_matrix(alg)
  # align pseudotime / weights to cells actually present in this TFA matrix
  cells <- colnames(P)

  pt_use <- pt_joint[cells, , drop=FALSE]
  cw_use <- cw_joint[cells, , drop=FALSE]

  keep_tf <- apply(P, 1, function(x) {all(is.finite(x)) && sd(x) > 1e-8})
  cat("Removing", sum(!keep_tf), "constant/invalid TFs\n")
  P <- P[keep_tf, , drop=FALSE]

  # fit continuous TFA with Gaussian GAM. Use nknots=6 for simplicity (otherwise takes too long)
  fit <- fitGAM(
    counts=P,
    pseudotime=pt_use,
    cellWeights=cw_use,
    nknots=6,
    family="gaussian",
    offset=rep(0, ncol(P)),
    verbose=TRUE
  )

  TFA_fits[[alg]] <- fit

  # within-lineage dynamics
  TFA_assoc[[alg]] <- as.data.frame(
    associationTest(
      fit,
      global=TRUE,
      lineages=TRUE
    )
  )

  # start vs end within each lineage
  TFA_start_end[[alg]] <- as.data.frame(
    startVsEndTest(
      fit,
      global=TRUE,
      lineages=TRUE
    )
  )

  # between-lineage trajectory-pattern differences
  TFA_pattern[[alg]] <- as.data.frame(
    patternTest(
      fit,
      global=TRUE,
      pairwise=TRUE
    )
  )

  # between-lineage endpoint differences
  TFA_end[[alg]] <- as.data.frame(
    diffEndTest(
      fit,
      global=TRUE,
      pairwise=TRUE
    )
  )
}

# save
save_tradeseq_results <- function(res, fit, P, outdir, type,
                                  lineage_names=colnames(pt_joint),
                                  top_n=10, plot_n=5) {

  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

  write.table(res, paste0(outdir, "/results.txt"),
              sep="\t", quote=FALSE, col.names=NA)

  # overall
  top10 <- data.frame(
    overall=rownames(res)[order(res$waldStat, decreasing=TRUE)[1:top_n]]
  )

  if (type == "lineage") {

    for (j in seq_along(lineage_names)) {

      # associationTest uses waldStat_1; startVsEndTest uses waldStat_lineage1
      candidates <- c(
        paste0("waldStat_", j),
        paste0("waldStat_lineage", j)
      )

      wcol <- candidates[candidates %in% colnames(res)][1]

      if (is.na(wcol))
        stop("Cannot find lineage Wald column for lineage ", j)

      o <- order(res[[wcol]], decreasing=TRUE)
      top10[[lineage_names[j]]] <- rownames(res)[o[1:top_n]]
    }

  } else if (type == "pairwise") {

    pairs <- combn(seq_along(lineage_names), 2)

    for (k in seq_len(ncol(pairs))) {
      i <- pairs[1,k]
      j <- pairs[2,k]

      wcol <- paste0("waldStat_", i, "vs", j)
      nm <- paste0(lineage_names[i], "_vs_", lineage_names[j])

      o <- order(res[[wcol]], decreasing=TRUE)
      top10[[nm]] <- rownames(res)[o[1:top_n]]
    }
  }

  write.table(
    top10,
    paste0(outdir, "/top10_genes.txt"),
    sep="\t", quote=FALSE, row.names=FALSE
  )

  # top 5 plots
  for (nm in names(top10)) {

    pdir <- paste0(outdir, "/", nm)
    dir.create(pdir, recursive=TRUE, showWarnings=FALSE)

    for (i in seq_len(min(plot_n, nrow(top10)))) {
      gene <- top10[[nm]][i]

      pdf(paste0(pdir, "/", i, "_", gene, ".pdf"), 7, 5)
      print(plotSmoothers(fit, P, gene=gene))
      dev.off()
    }
  }
}

dir.create("plots/tradeseq/TFA", recursive=TRUE, showWarnings=FALSE)

for (alg in nonnegative_algos) {

  cat("\nSaving:", alg, "\n")

  P <- get_TFA_matrix(alg)

  # same filtering used before fitting
  keep_tf <- apply(P, 1, function(x) all(is.finite(x)) && sd(x) > 1e-8)
  P <- P[keep_tf, , drop=FALSE]

  fit <- TFA_fits[[alg]]
  adir <- paste0("plots/tradeseq/TFA/", alg)

  save_tradeseq_results(
    TFA_assoc[[alg]], fit, P,
    paste0(adir, "/association"),
    type="lineage"
  )

  save_tradeseq_results(
    TFA_start_end[[alg]], fit, P,
    paste0(adir, "/startVsEnd"),
    type="lineage"
  )

  save_tradeseq_results(
    TFA_pattern[[alg]], fit, P,
    paste0(adir, "/pattern"),
    type="pairwise"
  )

  save_tradeseq_results(
    TFA_end[[alg]], fit, P,
    paste0(adir, "/diffEnd"),
    type="pairwise"
  )
}

# > sessionInfo()                                                                                                                                                     
# R version 4.3.3 (2024-02-29)                                                                                                                                        
# Platform: x86_64-conda-linux-gnu (64-bit)                                                                                                                           
# Running under: Rocky Linux 8.10 (Green Obsidian)                                                                                                                    
                                                                                                                                                                    
# Matrix products: default                                                                                                                                            
# BLAS/LAPACK: /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/benchmarking_tfa_algorithms/algorithms/BITFAM/lib/libopenblasp-r0.3.27.so;  LAPACK version
#  3.12.0   
# Random number generation:                                                                                                                                           
#  RNG:     Mersenne-Twister                                                                                                                                          
#  Normal:  Inversion                                                                                                                                                 
#  Sample:  Rounding                                                                                                                                                  
                                                                                                                                                                    
# locale:                                                                                                                                                             
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C                                                                                                                        
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8                                                                                                              
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8                                                                                                             
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                                                                                                                           
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C                                                                                                                      
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C                                                                                                                 
                                                                                                                                                                    
# time zone: America/Chicago                                                                                                                                          
# tzcode source: system (glibc)                                                                                                                                       
                                                                                                                                                                    
# attached base packages:                                                                                                                                             
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets                                                                                            
# [8] methods   base                                                                                                                                                  
                                                                                                                                                                    
# other attached packages:                                                                                                                                            
#  [1] Matrix_1.6-5                mgcv_1.9-1                                                                                                                         
#  [3] nlme_3.1-165                igraph_2.0.3                                                                                                                       
#  [5] Signac_1.14.0               RColorBrewer_1.1-3                                                                                                                 
#  [7] tradeSeq_1.16.0             presto_1.0.0                                                                                                                       
#  [9] data.table_1.15.4           Rcpp_1.0.13                                                                                                                        
# [11] OmnipathR_3.10.1            ggplot2_3.5.1                                                                                                                      
# [13] dplyr_1.1.4                 Seurat_5.1.0                                                                                                                       
# [15] SeuratObject_5.0.2          sp_2.1-4                                                                                                                           
# [17] slingshot_2.10.0            TrajectoryUtils_1.10.0                                                                                                             
# [19] SingleCellExperiment_1.24.0 SummarizedExperiment_1.32.0                                                                                                        
# [21] Biobase_2.62.0              GenomicRanges_1.54.1                                                                                                               
# [23] GenomeInfoDb_1.38.1         IRanges_2.36.0                                                                                                                     
# [25] S4Vectors_0.40.2            BiocGenerics_0.48.1                                                                                                                
# [27] MatrixGenerics_1.14.0       matrixStats_1.3.0                                                                                                                  
# [29] princurve_2.1.6                                                                                                                                                
                                                                                                                                                                    
# loaded via a namespace (and not attached):                                                                                                                          
#   [1] RcppAnnoy_0.0.22          splines_4.3.3                                                                                                                       
#   [3] later_1.3.2               bitops_1.0-9                                                                                                                        
#   [5] tibble_3.2.1              cellranger_1.1.0                                                                                                                    
#   [7] polyclip_1.10-7           fastDummies_1.7.3                                                                                                                   
#   [9] lifecycle_1.0.4           edgeR_3.42.4                                                                                                                        
#  [11] globals_0.16.3            lattice_0.22-6                                                                                                                      
#  [13] vroom_1.6.5               MASS_7.3-60.0.1                                                                                                                     
#  [15] backports_1.5.0           magrittr_2.0.3                                                                                                                      
#  [17] limma_3.58.1              plotly_4.11.0                                                                                                                       
#  [19] rmarkdown_2.27            yaml_2.3.10                                                                                                                         
#  [21] httpuv_1.6.15             sctransform_0.4.1                                                                                                                   
#  [23] spam_2.11-0               spatstat.sparse_3.1-0                                                                                                               
#  [25] reticulate_1.43.0         cowplot_1.1.3                                                                                                                       
#  [27] pbapply_1.7-2             lubridate_1.9.4                                                                                                                     
#  [29] abind_1.4-5               zlibbioc_1.48.0                                                                                                                     
#  [31] rvest_1.0.4               Rtsne_0.17   
#  [33] purrr_1.0.2               RCurl_1.98-1.16          
#  [35] rappdirs_0.3.3            GenomeInfoDbData_1.2.11  
#  [37] ggrepel_0.9.6             irlba_2.3.5.1            
#  [39] listenv_0.9.1             spatstat.utils_3.0-5     
#  [41] goftest_1.2-3             RSpectra_0.16-2          
#  [43] spatstat.random_3.3-1     fitdistrplus_1.2-1       
#  [45] parallelly_1.41.0         DelayedMatrixStats_1.24.0
#  [47] RcppRoll_0.3.1            leiden_0.4.3.1           
#  [49] codetools_0.2-20          DelayedArray_0.28.0      
#  [51] xml2_1.3.6                tidyselect_1.2.1         
#  [53] farver_2.1.2              viridis_0.6.5            
#  [55] spatstat.explore_3.3-1    jsonlite_1.8.8           
#  [57] progressr_0.15.1          ggridges_0.5.6           
#  [59] survival_3.8-3            systemfonts_1.1.0        
#  [61] tools_4.3.3               progress_1.2.3           
#  [63] ragg_1.3.3                ica_1.0-3                
#  [65] glue_1.7.0                gridExtra_2.3            
#  [67] SparseArray_1.2.2         xfun_0.46                
#  [69] withr_3.0.0               BiocManager_1.30.23      
#  [71] fastmap_1.2.0             fansi_1.0.6              
#  [73] digest_0.6.36             timechange_0.3.0         
#  [75] R6_2.5.1                  mime_0.12                
#  [77] textshaping_0.4.0         colorspace_2.1-1         
#  [79] scattermore_1.2           tensor_1.5               
#  [81] spatstat.data_3.1-2       utf8_1.2.4               
#  [83] tidyr_1.3.1               generics_0.1.3           
#  [85] prettyunits_1.2.0         httr_1.4.7               
#  [87] htmlwidgets_1.6.4         S4Arrays_1.2.0           
#  [89] uwot_0.2.3                pkgconfig_2.0.3          
#  [91] gtable_0.3.5              lmtest_0.9-40            
#  [93] XVector_0.42.0            selectr_0.4-2            
#  [95] htmltools_0.5.8.1         dotCall64_1.2            
#  [97] scales_1.3.0              png_0.1-8                
#  [99] spatstat.univar_3.0-0     knitr_1.48               
# [101] tzdb_0.4.0                reshape2_1.4.4           
# [103] checkmate_2.3.2           curl_5.2.1               
# [105] zoo_1.8-14                stringr_1.5.1            
# [107] KernSmooth_2.23-26        miniUI_0.1.1.1           
# [109] pillar_1.9.0              grid_4.3.3               
# [111] logger_0.4.0              vctrs_0.6.5              
# [113] RANN_2.6.2                promises_1.3.0           
# [115] xtable_1.8-4              cluster_2.1.8            
# [117] evaluate_0.24.0           readr_2.1.5              
# [119] Rsamtools_2.18.0          locfit_1.5-9.10          
# [121] cli_3.6.3                 compiler_4.3.3           
# [123] rlang_1.1.4               crayon_1.5.3             
# [125] future.apply_1.11.3       labeling_0.4.3           
# [127] plyr_1.8.9                stringi_1.8.4            
# [129] BiocParallel_1.36.0       viridisLite_0.4.2        
# [131] deldir_2.0-4              Biostrings_2.70.1        
# [133] munsell_0.5.1             lazyeval_0.2.2           
# [135] spatstat.geom_3.3-2       RcppHNSW_0.6.0           
# [137] hms_1.1.3                 patchwork_1.3.0          
# [139] sparseMatrixStats_1.14.0  bit64_4.5.2              
# [141] future_1.34.0             statmod_1.5.0            
# [143] shiny_1.8.1.1             ROCR_1.0-11              
# [145] fastmatch_1.1-6           bit_4.5.0.1              
# [147] readxl_1.4.3   