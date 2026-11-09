##
#### This script is adapted from
## ---------------------------
## Project name: TIGER
## Purpose of script: Fig 4 & Fig S4-7
## Author: Chen Chen
## Date Created: 2022-10-28
##
## Copyright (c) Chen Chen, 2022
## Email: cchen22@arizona.edu

##### Part I. Make prior, expression, gene names, tree, and tf names
library(Seurat)
# library(Signac)
# library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(OmnipathR)
library(slingshot)
library(parallel)

max_cores = 50

set.seed(1234) # 

# Create Seurat object. from TIGER's zenodo dir
pbmc <- readRDS("TIGER_zenodo/data/PBMC/pbmc_final.rds")

# load prior
collectri_grn <- collectri(organism=9606L, genesymbols=TRUE, loops=TRUE)


# library(JASPAR2020)
# library(TFBSTools)
# library(motifmatchr)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(presto)

markers_rna <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(pbmc, id = motif.names)


## filter expression by marker genes and prior by marker motifs
padj.cutoff = 0.01
logFC.cutoff = 0

ctmarkers_rna <- markers_rna %>%
  dplyr::filter(RNA.padj < padj.cutoff, RNA.logFC > logFC.cutoff) %>% 
  dplyr::arrange(-RNA.auc) %>%
  dplyr::select(gene,RNA.group)

ctmarkers_motif <- markers_motifs[markers_motifs$motif.padj < 1e-2, ]
ctmarkers_motif <- ctmarkers_motif[ctmarkers_motif$gene %in% collectri_grn$source_genesymbol, ] # filter regulators without prior out edges (ambiguity in column permutation)
ctmarkers_motif <- ctmarkers_motif[ctmarkers_motif$gene %in% rownames(pbmc$SCT), ] # must be expressed in some cells

## export expression, normalized expression, prior, and cell names
normalized_expr = GetAssayData(pbmc,assay = "SCT",slot = "data") 
expr = GetAssayData(pbmc,assay = "RNA",slot = "data") 
filtered_normalized_expr <- normalized_expr[rownames(normalized_expr) %in% ctmarkers_rna$gene, ]
filtered_expr <- expr[rownames(expr) %in% ctmarkers_rna$gene, ]

system('mkdir -p TIGER_zenodo/data/PBMC_ct')


cell_counts_per_c <- mclapply(unique(pbmc$celltype), function(c){
  cur_expr <- filtered_expr[,colnames(filtered_expr) %in% names(pbmc$celltype[pbmc$celltype==c])]
  write.table(cur_expr, paste0('TIGER_zenodo/data/PBMC_ct/', gsub(' ', '_', c), '_exp.txt'), sep='\t', quote=FALSE, col.names = FALSE)
  write.table(colnames(cur_expr), paste0('TIGER_zenodo/data/PBMC_ct/', gsub(' ', '_', c), '_cellnames.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
  cur_expr <- filtered_normalized_expr[,colnames(filtered_normalized_expr) %in% names(pbmc$celltype[pbmc$celltype==c])]
  write.table(cur_expr, paste0('TIGER_zenodo/data/PBMC_ct/normalized_', gsub(' ', '_', c), '_exp.txt'), sep='\t', quote=FALSE, col.names = FALSE)
  return(ncol(cur_expr))
  }, mc.cores=max_cores)

cell_counts_per_c <- unlist(cell_counts_per_c)
names(cell_counts_per_c) <- gsub(" ", "_", unique(pbmc$celltype))

filtered_collectri_grn <- collectri_grn[collectri_grn$source_genesymbol %in% ctmarkers_motif$gene & collectri_grn$target_genesymbol %in% rownames(normalized_expr), ]

filtered_collectri_grn <- filtered_collectri_grn %>%
  select(source_genesymbol, target_genesymbol) %>%  
  mutate(conf = 1) 

write.table(filtered_collectri_grn, "TIGER_zenodo/data/PBMC_ct/pbmc10k_prior.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(filtered_collectri_grn, "TIGER_zenodo/pbmc10k_prior.txt", sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

system('mkdir -p TIGER_zenodo/data/PBMC_ct/all_exps')

write.table(rownames(filtered_expr), paste0('TIGER_zenodo/data/PBMC_ct/all_genes.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
write.table(unique(ctmarkers_motif$gene), paste0('TIGER_zenodo/data/PBMC_ct/all_regulators.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)

workdir <- system('pwd', intern=TRUE)

all_priors <- data.frame(celltype=gsub(" ", "_", unique(pbmc$celltype)))
all_priors$path <- paste0(workdir, '/TIGER_zenodo/pbmc10k_prior.txt')
all_priors$conf <- 1
write.table(all_priors, "TIGER_zenodo/data/PBMC_ct/all_exps/all_priors.txt", sep = '\t', quote=FALSE, col.names=FALSE, row.names=FALSE)

all_exp <- data.frame(celltype=gsub(" ", "_", unique(pbmc$celltype)))
all_exp$path <- paste0(workdir, '/TIGER_zenodo/data/PBMC_ct/normalized_', all_exp$celltype, '_exp.txt')
write.table(all_exp, "TIGER_zenodo/data/PBMC_ct/all_exps/all_exps.txt", sep = '\t', quote=FALSE, col.names=FALSE, row.names=FALSE)

# make cell lineage tree from the weighted nearest neighbor embedding
dimred <- pbmc@reductions$wnn.umap@cell.embeddings
clustering <- pbmc$celltype
lineages <- getLineages(data = dimred, clusterLabels = clustering, start.clus="HSPC")

# Plot the lineages
png("lineages_plot.png", width = 3200, height = 1600, res = 300) 
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
    text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(SlingshotDataSet(lineages), lwd = 3, col = "black")
dev.off()

lineage_dat <- as.list(slingLineages(lineages))


lineage_tree <- list()
for(l in lineage_dat){
  for(i in 1:(length(l) - 1)){
    lineage_tree[[length(lineage_tree)+1]] <- c(l[[i+1]], l[[i]])
  }
}
lineage_tree <- as.data.frame(gsub(" ", "_", do.call(rbind, unique(lineage_tree))))
lineage_tree[, 3] <- 1
write.table(lineage_tree, "TIGER_zenodo/data/PBMC_ct/all_exps/celltype_tree_ancestor.txt", sep = '\t', quote=FALSE, col.names=FALSE, row.names=FALSE)




# ##### Part II. After running all algorithms, compare the results with chromVAR deviation as a genome-occupancy based TFA proxy
chrom_assay = GetAssayData(pbmc,assay = "chromvar",slot = "data") 
filtered_chrom_assay <- chrom_assay[rownames(chrom_assay) %in% ctmarkers_motif$motif.feature,]

# we can do this for the special case that motif and TF are of 1 to 1 correspondence
rownames(filtered_chrom_assay) <- sapply(rownames(filtered_chrom_assay), function(x){(ctmarkers_motif[ctmarkers_motif$motif.feature==x,])[['gene']][[1]]})
write.table(filtered_chrom_assay, "chromVAR_TFA.txt", quote=FALSE, row.names=TRUE, col.names=TRUE, sep = '\t')



# ##### Deprecated Part III. plots start here ######

# ## Fig 4C. Heatmap of top5 TFs with all cells
# heat.data = readRDS(paste0("result/PBMC/TIGER_mean.rds"))
# tiger_top5 = readRDS("result/PBMC/tiger_top5.rds")

# top_TFs = unique(unlist(tiger_top5))
# heat.data = apply(heat.data,2,function(x) x/max(x))
# heat.data = heat.data[top_TFs,]

# library(pheatmap)
# library(ggsci)
# library(RColorBrewer)
# library(scales)
# annotate_col = data.frame(CellType = as.factor(colnames(heat.data)))
# rownames(annotate_col) = colnames(heat.data)
# ann_colors = list(show_col(hue_pal()(19)))

# heat = pheatmap(
#   heat.data,
#   scale = 'none',
#   color = brewer.pal(9,"Oranges"),
#   cluster_rows = F,
#   cluster_cols = F,
#   show_colnames = F,
#   show_rownames = T,
#   treeheight_row = 0,
#   treeheight_col = 0,
#   border_color = 'grey',
#   annotation_col = annotate_col,
#   annotation_colors = ann_colors
# )
# ggsave("figure/fig_4C.png",heat,width = 6,height = 8)

# ## Fig 4. UMAP
# pbmc = readRDS("data/PBMC/pbmc_final.rds")

# ## cell type annotation: Fig 4A
# pbmc@meta.data$celltype = factor(pbmc@meta.data$celltype,levels = colnames(heat.data))
# celltype = DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", 
#                    label = TRUE, label.size = 4, repel = TRUE) + 
#   ggtitle("Cell Type")+
#   theme(legend.position = "none")
# ggsave("figure/fig_4A.png",celltype,width = 6,height = 4)

# ## RNA expression feature plot: Fig 4D-F, Fig S4-6
# IRF4_RNA = FeaturePlot(pbmc, features = "sct_IRF4", reduction = 'wnn.umap',
#                        label = T,label.size = 2.5,repel = T) + ggtitle("RNA_IRF4")
# ggsave("figure/IRF4_RNA.png",IRF4_RNA,width = 4,height = 4)

# TBX21_RNA = FeaturePlot(pbmc, features = "sct_TBX21", reduction = 'wnn.umap',
#                         label = T,label.size = 2.5,repel = T) + ggtitle("RNA_TBX21")
# ggsave("figure/TBX21_RNA.png",TBX21_RNA,width = 4,height = 4)

# TCF4_RNA = FeaturePlot(pbmc, features = "sct_TCF4", reduction = 'wnn.umap',
#                        label = T,label.size = 2.5,repel = T) + ggtitle("RNA_TCF4")
# ggsave("figure/TCF4_RNA.png",TCF4_RNA,width = 4,height = 4)

# PAX5_RNA = FeaturePlot(pbmc, features = "sct_PAX5", reduction = 'wnn.umap',
#                        label = T,label.size = 2.5,repel = T) + ggtitle("RNA_PAX5")
# ggsave("figure/PAX5_RNA.png",PAX5_RNA,width = 4,height = 4)

# ## Multi-modal Feature plot: Fig 4D-F, Fig S4-6
# DefaultAssay(pbmc) = "ATAC"

# motif.name <- ConvertMotifID(pbmc, name = 'IRF4')
# IRF4_motif = FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, reduction = 'wnn.umap',
#                          label = T,label.size = 2.5,repel = T) + ggtitle("motif_IRF4")
# ggsave("figure/IRF4_motif.png",IRF4_motif,width = 4,height = 4)

# motif.name <- ConvertMotifID(pbmc, name = 'TBX21')
# TBX21_motif = FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, reduction = 'wnn.umap',
#                           label = T,label.size = 2.5,repel = T) + ggtitle("motif_TBX21")
# ggsave("figure/TBX21_motif.png",TBX21_motif,width = 4,height = 4)

# motif.name <- ConvertMotifID(pbmc, name = 'TCF4')
# TCF4_motif = FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, reduction = 'wnn.umap',
#                          label = T,label.size = 2.5,repel = T) + ggtitle("motif_TCF4")
# ggsave("figure/TCF4_motif.png",TCF4_motif,width = 4,height = 4)

# motif.name <- ConvertMotifID(pbmc, name = 'PAX5')
# PAX5_motif = FeaturePlot(pbmc, features = motif.name, min.cutoff = 0, reduction = 'wnn.umap',
#                          label = T,label.size = 2.5,repel = T) + ggtitle("motif_PAX5")
# ggsave("figure/PAX5_motif.png",PAX5_motif,width = 4,height = 4)

# ## TFA feature plots: Fig 4D-F, Fig S4-6
# DefaultAssay(pbmc) = "TIGER"
# IRF4_TF = FeaturePlot(pbmc, features = "IRF4", reduction = 'wnn.umap',
#                       label = T,label.size = 2.5,repel = T)+ ggtitle("TF_IRF4")
# ggsave("figure/IRF4_TF.png",IRF4_TF,width = 4,height = 4)

# TBX21_TF = FeaturePlot(pbmc, features = "TBX21", reduction = 'wnn.umap',
#                        label = T,label.size = 2.5,repel = T) + ggtitle("TF_TBX21")
# ggsave("figure/TBX21_TF.png",TBX21_TF,width = 4,height = 4)

# TCF4_TF = FeaturePlot(pbmc, features = "TCF4", reduction = 'wnn.umap',
#                       label = T,label.size = 2.5,repel = T) + ggtitle("TF_TCF4")
# ggsave("figure/TCF4_TF.png",TCF4_TF,width = 4,height = 4)

# PAX5_TF = FeaturePlot(pbmc, features = "PAX5", reduction = 'wnn.umap',
#                       label = T,label.size = 2.5,repel = T) + ggtitle("TF_PAX5")
# ggsave("figure/PAX5_TF.png",PAX5_TF,width = 4,height = 4)




# ## Fig 4B. comparing TIGER vs. VIPER consistency with scATACseq
# # RBO
# library(gespeR)
# library(ggpubr)
# tiger_rank = readRDS("result/PBMC/tiger_rank.rds")
# viper_rank = readRDS("result/PBMC/viper_rank.rds")
# motif_rank = readRDS("result/PBMC/motif_rank.rds")

# rbo_tg = mapply(gespeR::rbo,motif_rank,tiger_rank,p=1,side='top')
# rbo_vp = mapply(gespeR::rbo,motif_rank,viper_rank,p=1,side='top')

# x = rep(names(rbo_tg),2)
# y = c(rbo_tg,rbo_vp)
# group = c(rep("TIGER",19),rep("VIPER",19))
# dfm = data.frame('x'=x,'y'=y,'group'=group)

# dot.plot = ggdotchart(dfm, x = "x", y = "y",
#                       color = "group",                               
#                       palette = 'npg', # Custom color palette
#                       sorting = "none",                      
#                       rotate = F,                               
#                       dot.size = 4,                                 
#                       y.text.col = F,   
#                       position = position_dodge(0.8),
#                       alpha = 0.2,
#                       add = "segment",
#                       ylab = "RBO similarity",
#                       xlab = "",
#                       ggtheme = theme_pubr()                       
# )+
#   #theme_cleveland()+
#   theme(legend.title = element_blank())+
#   geom_hline(yintercept = mean(rbo_tg), linetype = "dashed", color = "#E64B35FF" )+
#   geom_hline(yintercept = mean(rbo_vp), linetype = "dashed", color = "#4DBBD5FF" )

# ggsave("figure/fig_4B.png",dot.plot,width = 6,height = 4)

# wilcox.test(rbo_tg,rbo_vp,alternative = "greater") # p = 0.03179


# ## Figure S7. edge consitency with ChIPseq
# TF_chip = read.table(paste0("data/PBMC/cistrome/","IRF4_Bcell",".txt"))
# colnames(TF_chip) = c("chrom",	"txStart",	"txEnd",	"refseq",	"score",	"strand",	"symbol")
# TF_chip = TF_chip %>%
#   group_by(symbol)%>%
#   summarise("max.score" = median(score)) %>%
#   arrange(desc(max.score)) %>%
#   pull(max.score,symbol)

# ss = readRDS(paste0("result/PBMC/hpc/TIGER_","Intermediate B",".rds"))

# TF_tiger = ss$W[,'IRF4']
# TF_tiger = TF_tiger[TF_tiger!=0]

# TGs = intersect(names(TF_tiger),names(TF_chip))
# TF_chip = TF_chip[TGs]
# TF_tiger = abs(TF_tiger[TGs])

# # plot
# edge_df = data.frame(x=TF_tiger,y=TF_chip,name=names(TF_tiger))
# edge_cor = ggscatter(edge_df, x = "x", y = "y",
#                      xlab = "TIGER edge weight",
#                      ylab = "ChIP edge score",
#                      cor.coef = T,
#                      cor.coeff.args = list(method = "kendall", label.x.npc = "middle", label.y.npc = "top"),
#                      color = "black", 
#                      size = 1, alpha = 0.6,
#                      position=position_jitter(h=0.8,w=0.8),
#                      add = "reg.line",
#                      #label = "name", 
#                      repel = TRUE,
#                      #label.select = list(criteria = "`x` > 10 & `y` < 10")
# )
# ggsave("plot/PBMC10K/edge_cor.png",edge_cor,width = 4,height = 4)


# bar_df = data.frame(x = names(TF_tiger),y=TF_tiger)
# bar.plot = ggbarplot(bar_df, x="x", y="y",
#                      color = "white",
#                      fill = 'red',
#                      palette ='npg',
#                      sort.val = "asc",
#                      alpha = 0.5,
#                      position = position_dodge2(padding = 0),
#                      xlab = "",
#                      ylab = "TIGER edge weight of IRF4",
#                      legend.title="",
# )+coord_flip()
# ggsave("plot/PBMC10K/bar_plot.png",bar.plot,width = 4,height = 4)


