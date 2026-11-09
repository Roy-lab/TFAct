## Adapted from:
## /mnt/dv/wid/projects7/Roy-singlecell2/sr_work/netinf/scMTNI_project/Rscripts/makeTopicRegBubble_ggplot_SZ_FBS.R
##

args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(data.table)
library(tidyr)
library(dplyr)
library(gridExtra)

# my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")        # yellow-green-blue

#PREFIX FOR PLOTS
# codename_prefix="/mnt/dv/wid/projects5/Roy-singlecell/saptarshi_work/scMTNImanuscript/results/scMTNI/a2s/condorOutputs/analysis/lda_TFcellbygene/network_cf0.8/k10/downstream_analysis/"
codename_prefix=args[1] # /mnt/dv/wid/projects6/Roy-singlecell3/sridharanlab2/Andrew_multiome/spencer_work/results/condor_results/sparse_-4_prior_4/analysis/lda_TFcellbygene/k10/
num_topics=args[2] # 10
xordfile <- args[3] # /mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/dynamic_networks/Andrew_multiome_project/full_tree/cluster_id.txt
dprop_list <- list()


datatab1=NULL
# datatab1_TFA=NULL
# datatab1_mRNA=NULL

# rtop=5

## For regulator nodes with out-degree > "maxdeg,"
## reduce their out-degree to "maxdeg."
## It is done to upper bound the size of the
## bubbles corr. to these nodes so that
## the bubble plot looks more aesthetic.
maxdeg=40

## there may be too many regulators in a topic. 
## get top regulators by the sum of ntarget over 
## cell clusters. We include all tied regulators
max_regs_in_a_plot = 10

# # the following specifies the minimum out degree, summed over all cell clusters
# mRNA_filter=5
# TFA_filter=7 # at least 7 total targets across all cell types are needed to be plotted in the TFA specific bubble plot
overall_filter=10


## (deprecated) the following further filters regulators such that they must have at least x number of targets in at least 1 cell cluster
# mRNA_max_c_min_ntarget = 3
# TFA_max_c_min_ntarget = 5

if (length(args) > 3) {
  outpath=args[4]
} else {
  outpath=codename_prefix
}

if(length(args) > 4){
  indir=args[5] #"/mnt/dv/wid/projects7/Roy-singlecell2/suvo_work/leukemic_patient_specific_study_updated_Oct2025/Patient5_PB3/scMTNI/data"
  cells <- read.csv(xordfile,header=F, sep = "\t")$V1
  for (cell in cells) {
    dt <- fread(paste0(indir, '/', cell, ".table"))
    dt$Gene <- gsub(paste0("_", cell), "", dt$Gene)
    
    ncell_withexpr <- rowSums(dt[, -1] > 0)
    prop_withexpr <- ncell_withexpr / (ncol(dt) - 1)
    
    dprop_list[[cell]] <- data.frame(
      regulator = dt$Gene,
      prop = prop_withexpr
    ) %>%
      data.table::setnames("prop", cell)
  }
  dprop_withexpr <- Reduce(function(x, y) merge(x, y, by="regulator", all=TRUE), dprop_list)

  # replace NAs with 0
  dprop_withexpr[is.na(dprop_withexpr)] <- 0

  dprop_withexpr_long=dprop_withexpr  %>% pivot_longer(
    cols = 2:ncol(dprop_withexpr),
    names_to = c("cell"),
    values_to = "prop"
  )
}



## For each topic
for(i in 1:as.numeric(num_topics)){
  datafile <- paste(c(codename_prefix,"/networks_lda_k", num_topics,"_topic", i, ".txt"),collapse="")

  # xordfile <- paste(c(codename_prefix,"cellorder.txt"),collapse="")
  ## Use shortened version of the cluster IDs, e.g., use "C8" for "cluster 8."
  ## However, the output is gibberish when short IDs are used. So, use the long names only.
  # xordfile <- "/mnt/dv/wid/projects5/Roy-singlecell/saptarshi_work/scMTNImanuscript/data/sridharan_a2s/celltype_order_short.txt"
  # yordfile <- paste(c(codename_prefix,"topic",i,"_min10_yorder.txt"),collapse="")

  data <- read.csv(datafile,header=T,sep='\t')
  data <- data[data$regulator != data$target,] # remove self loops
  datatab <- as.data.table(data[,c('regulator', 'cell')])[, .N, by= c('regulator', 'cell')][order(regulator)]
  colnames(datatab) <- c('yaxis', 'xaxis', 'ntarget')
  ## For regulator nodes with out-degree > "maxdeg,"
  ## reduce their out-degree to "maxdeg."
  ## It is done to upper bound the size of the 
  ## bubbles corr. to these nodes so that
  ## the bubble plot looks more aesthetic.
  ids=which(datatab$ntarget>maxdeg)
  # datatab$Degree[ids]=maxdeg;
  datatab$ntarget[ids]=maxdeg;

  xord <- read.csv(xordfile,header=F, sep = "\t")
  datatab$xaxis <- factor(datatab$xaxis,levels=xord[[1]])

  ## the following orders regulators the same way as /mnt/dv/wid/projects5/Roy-singlecell/saptarshi_work/scMTNImanuscript/code/scMTNI/Scripts/Network_Analysis/componentsanalysis_directed/getDirectedComponents
  ## by the number of unique targets (regardless of cell types).
  yord_target_counts <- sort(table(unique(data[,c('regulator', 'target')])$regulator)) # one may filter based on this

  yord <- names(yord_target_counts)
  datatab$yaxis <- factor(datatab$yaxis,levels=yord)

  datatab$topic=i
  datatab1=rbind(datatab1,datatab)

  # datatab_TFA = datatab[grepl("_TFA$", datatab$yaxis),]
  # datatab_mRNA = datatab[!grepl("_TFA$", datatab$yaxis),]

  # datatab_mRNA = datatab_mRNA[yaxis %in% datatab_mRNA[, .(val = sum(ntarget)), by = yaxis][val>=mRNA_filter,][, yaxis],]
  # datatab_mRNA = datatab_mRNA[yaxis %in% datatab_mRNA[datatab_mRNA$ntarget>=mRNA_max_c_min_ntarget,]$yaxis,]

  # datatab_TFA = datatab_TFA[yaxis %in% datatab_TFA[, .(val = sum(ntarget)), by = yaxis][val>=TFA_filter,][, yaxis],]
  # datatab_TFA = datatab_TFA[yaxis %in% datatab_TFA[datatab_TFA$ntarget>=TFA_max_c_min_ntarget,]$yaxis,]

  # datatab1_mRNA=rbind(datatab1_mRNA,datatab_mRNA)
  # datatab1_TFA=rbind(datatab1_TFA,datatab_TFA)
}


##############################################################
## Begin: Generate a bubble plot suitable for the manuscript
###############################################################
# # outname=paste0(outpath,"alltopics_rtop",rtop,"maxdeg",maxdeg,"_regulators_min10_bubbleplot_v2.pdf")

outname=paste0(outpath,"/alltopics_", "maxreg", max_regs_in_a_plot, "_maxdeg",maxdeg,"_regulators_min", overall_filter, "_bubbleplot.pdf")
dynamic_height = as.numeric(num_topics) * 0.8
if (max_regs_in_a_plot < 15){
  dynamic_height = as.numeric(num_topics) * 0.6
}
cells_plot <- sub("^cluster([0-9]+)$", "C\\1", cells)
pdf(outname, height = dynamic_height, width = 15,useDingbats=FALSE)

# ## yellow-green-blue
my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")
datatab1_overall = datatab1[yaxis %in% datatab1[, .(val = sum(ntarget)), by = yaxis][val>=overall_filter,][, yaxis],]

datatab1_overall_filtered <- datatab1_overall[, {
      yval <- .SD[, .(val = sum(ntarget)), by = yaxis]
      if (nrow(yval) <= max_regs_in_a_plot) {
          .SD
      } else {
          cutoff_val <- sort(yval$val, decreasing = TRUE)[max_regs_in_a_plot]
          keep_y <- yval[val >= cutoff_val, yaxis]
          .SD[yaxis %in% keep_y]
      }
    },by = topic]

setDT(datatab1_overall_filtered)
if(length(args) > 4){
  setDT(dprop_withexpr_long)

  # Regulators selected from edge-based filtering
  regulators_for_prop <- unique(
    datatab1_overall_filtered[, .(topic, yaxis)]
  )

  # All topic-regulator x cell combinations
  prop_plot_tab <- regulators_for_prop[
    , .(xaxis = cells),
    by = .(topic, yaxis)
  ]

  # Matching key: X_TFA should use expression of X
  prop_plot_tab[, regulator_match := sub("_TFA$", "", yaxis)]

  # Attach expression proportion for every cell
  prop_plot_tab[
    dprop_withexpr_long,
    prop := i.prop,
    on = .(regulator_match = regulator, xaxis = cell)
  ]

  # Preserve cell order
  prop_plot_tab[, xaxis := sub("^cluster([0-9]+)$", "C\\1", xaxis)][, xaxis := factor(xaxis, levels = cells_plot)]
}

# Preserve desired cluster order
datatab1_overall_filtered[, xaxis := sub("^cluster([0-9]+)$", "C\\1", xaxis)][, xaxis := factor(xaxis, levels = cells_plot)]
# cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = Degree,color=Degree)) +
cdat_sp <- ggplot(datatab1_overall_filtered, aes(x = xaxis, y = yaxis, size = ntarget,color = ntarget)) +
  geom_point()+
  # ggtitle(label = "top regulators")+
  ggtitle(label = paste0("Regulators with out degree >= ", overall_filter, "\nShowing top ", max_regs_in_a_plot, " regulators per topic"))+
  xlab(element_blank())+
  ylab(element_blank())+
  # theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$Degree)))+
  theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$ntarget)))+
  # theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
  theme(legend.position="bottom")+
  # theme(text = element_text(size=8),axis.text = element_text(size = 8),axis.text.x = element_text(size = 6.5),
  theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 90),
        strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=5, scales="free_y",labeller = "label_both")

plot(cdat_sp);
dev.off()

### cell proportion expressed for each regulator
if(length(args) > 4){
  outname <- paste0(
  outpath,
  "/alltopics_maxreg", max_regs_in_a_plot,
  "_maxdeg", maxdeg,
  "_regulators_min", overall_filter,
  "_bubbleplot_prop.pdf"
)

pdf(outname, height = dynamic_height, width = 15, useDingbats = FALSE)

my_cols <- c("#2166AC", "#B2182B")



cdat_sp <- ggplot(
  prop_plot_tab,
  aes(x = xaxis, y = yaxis, size = prop, color = prop)
) +
  geom_point() +
  ggtitle(
    "Proportion of cells in each cell cluster with non-zero expression for top regulators"
  ) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  scale_color_gradientn(
    colors = my_cols,
    limits = c(0, 1)
  ) +
  scale_size(range = c(0.1, 5)) +
  facet_wrap(
    ~topic,
    ncol = 5,
    scales = "free_y",
    labeller = "label_both"
  ) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 8),
    axis.text.y = element_text(size = 6.5),
    axis.text.x = element_text(size = 6.5, angle = 90),
    strip.text = element_text(size = 8)
  )

plot(cdat_sp)
dev.off()
}

############################################################
## End: Generate a bubble plot suitable for the manuscript
###########################################################

# for comparing against known list
allregs=paste0(outpath,"/alltopics_maxreg",max_regs_in_a_plot,"_maxdeg",maxdeg,"_regulators_min", overall_filter,  "_allregs.txt")
write.table(unique(datatab1_overall_filtered$yaxis), allregs, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(datatab1_overall_filtered, paste0(outpath,"/alltopics_maxreg",max_regs_in_a_plot,"_maxdeg",maxdeg,"_regulators_min", overall_filter, "_bubbleplot.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')


############################################################
## Begin: For _TFA regulators, identify other regulators sharing similar motifs
###########################################################

if(length(args) > 5){
motifIndir = args[6]
library(igraph)

## regulator set from bubble plot
regset <- unique(sub(
  "_TFA$", "",
  datatab1_overall_filtered[grepl("_TFA$", yaxis), yaxis]
))

## load families across clusters
fam_all <- rbindlist(lapply(cells, function(cell) {
  f <- fread(file.path(
    motifIndir,
    paste0(cell, "_motif_families_filteredexpression.txt")
  ))
  f[, .(TFs)]
}))

## keep families involving at least one displayed regulator
fam_all[, members := strsplit(TFs, ",\\s*")]
fam_all <- fam_all[
  sapply(members, function(z) any(z %in% regset))
]

## merge families whenever they overlap
edges <- unique(rbindlist(lapply(fam_all$members, function(z) {
  if (length(z) < 2) return(NULL)
  x <- t(combn(z, 2))
  data.table(a = x[,1], b = x[,2])
})))

g <- graph_from_data_frame(edges, directed = FALSE)
m <- components(g)$membership

families <- split(names(m), m)
families <- families[sapply(families, function(z) any(z %in% regset))]

## sort: displayed regulators first
families <- lapply(families, function(z)
  c(sort(intersect(z, regset)),
    sort(setdiff(z, regset)))
)

outname <- paste0(
  outpath, "/alltopics_",
  "maxreg", max_regs_in_a_plot,
  "_maxdeg", maxdeg,
  "_regulators_min", overall_filter,
  "_bubbleplot_motif_families.pdf"
)

names(families) <- paste0("family ", seq_along(families))

pdf(outname, width = 12, height = 11, family = "Helvetica")
par(mar = c(0.3, 0.3, 0.3, 0.3))

plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

cex <- 0.68
line_h <- 0.014
family_gap <- 0.004

xcols <- c(0.02, 0.51)
xmax  <- c(0.48, 0.97)

## approximately split families equally between columns
split_at <- ceiling(length(families) / 2)
cols <- list(
  seq_len(split_at),
  if (split_at < length(families))
    (split_at + 1):length(families)
  else integer(0)
)

for (cc in 1:2) {

  x0 <- xcols[cc]
  y <- 0.985

  for (i in cols[[cc]]) {

    fam <- families[[i]]

    ## family label
    label <- paste0(names(families)[i], ": ")
    text(x0, y, label, adj = c(0, 0.5), cex = cex, font = 2)

    x <- x0 + strwidth(label, cex = cex, font = 2)

    for (j in seq_along(fam)) {

      tf <- fam[j]
      shown <- tf %in% regset

      lab <- if (shown) paste0(tf, "_TFA") else tf
      if (j < length(fam))
        lab <- paste0(lab, ", ")

      w <- strwidth(lab, cex = cex)

      ## wrap within column
      if (x + w > xmax[cc]) {
        x <- x0
        y <- y - line_h
      }

      text(
        x, y, lab,
        adj = c(0, 0.5),
        cex = cex,
        col = if (shown) "blue" else "black"
      )

      x <- x + w
    }

    y <- y - line_h - family_gap
  }
}

dev.off()
}




##############################################################
## Begin: Generate TFA and mRNA specific plots
###############################################################
# outname=paste0(outpath,"alltopics_mRNA_","maxdeg",maxdeg,"_regulators_min", mRNA_filter,  "_", mRNA_max_c_min_ntarget,"_bubbleplot.pdf")
# pdf(outname, height = 8,width = 15,useDingbats=FALSE)

# ## yellow-green-blue
# my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")

# # cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = Degree,color=Degree)) +
# cdat_sp <- ggplot(datatab1_mRNA, aes(x = xaxis, y = yaxis, size = ntarget,color = ntarget)) +
#   geom_point()+
#   # ggtitle(label = "top regulators")+
#   ggtitle(label =  paste0("Regulators with total out degrees >= ", mRNA_filter, "\n","and have >= ", mRNA_max_c_min_ntarget, " targets in a cluster"))+
#   xlab(element_blank())+
#   ylab(element_blank())+
#   # theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$Degree)))+
#   theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$ntarget)))+
#   # theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
#   theme(legend.position="bottom")+
#   # theme(text = element_text(size=8),axis.text = element_text(size = 8),axis.text.x = element_text(size = 6.5),
#   theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 90),
#         strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=5, scales="free_y",labeller = "label_both")

# plot(cdat_sp);
# dev.off()

# outname=paste0(outpath,"alltopics_TFA_","maxdeg",maxdeg,"_regulators_min", TFA_filter, "_", TFA_max_c_min_ntarget, "_bubbleplot.pdf")
# pdf(outname, height = 8,width = 15,useDingbats=FALSE)

# ## yellow-green-blue
# my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")

# # cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = Degree,color=Degree)) +
# cdat_sp <- ggplot(datatab1_TFA, aes(x = xaxis, y = yaxis, size = ntarget,color = ntarget)) +
#   geom_point()+
#   # ggtitle(label = "top regulators")+
#   ggtitle(label = paste0("Regulators with total out degrees >= ", TFA_filter, "\n","and have >= ", TFA_max_c_min_ntarget, " targets in a cluster"))+
#   xlab(element_blank())+
#   ylab(element_blank())+
#   # theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$Degree)))+
#   theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$ntarget)))+
#   # theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
#   theme(legend.position="bottom")+
#   # theme(text = element_text(size=8),axis.text = element_text(size = 8),axis.text.x = element_text(size = 6.5),
#   theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 90),
#         strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=5, scales="free_y",labeller = "label_both")

# plot(cdat_sp);
# dev.off()
##############################################################
## End: Generate TFA and mRNA specific plots
###############################################################


#####################################################
## Begin: Generate a bubble plot suitable for slides
#####################################################
# outname=paste0(outpath,"alltopics_rtop",rtop,"maxdeg",maxdeg,"_regulators_min10_bubbleplot_v2.pdf")
# outname=paste0(outpath,"alltopics_","maxdeg",maxdeg,"_regulators_min10_bubbleplot_v2.pdf")
# # pdf(outname, height = 6,width = 15,useDingbats=FALSE)
# pdf(outname, height = 8,width = 15,useDingbats=FALSE)

# ## yellow-green-blue
# my_cols <-c("#001B87FF","#31DB92FF","#FCC201FF")

# # cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = Degree,color=Degree)) +
# cdat_sp <- ggplot(datatab1, aes(x = xaxis, y = yaxis, size = ntarget,color = ntarget)) +
#   geom_point()+
#   # ggtitle(label = "top regulators")+
#   ggtitle(label = "Regulators with out degree >= 10")+
#   xlab(element_blank())+
#   ylab(element_blank())+
#   # theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$Degree)))+
#   theme_bw()+scale_color_gradientn(colors = my_cols,limits=c(0,max(datatab1$ntarget)))+
#   # theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
#   # theme(legend.position="bottom")+
#   theme(legend.position="none")+
#   # theme(text = element_text(size=8),axis.text = element_text(size = 8),axis.text.x = element_text(size = 6.5),
#   # theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 90),
#   theme(text = element_text(size=8),axis.text.y = element_text(size = 6.5),axis.text.x = element_text(size = 6.5, angle = 15),
#         # strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=5, scales="free_y",labeller = "label_both")
#         strip.text = element_text(size = 8))+scale_size(range = c(0.1,5))+facet_wrap(~topic,ncol=4, scales="free_y",labeller = "label_both")

# plot(cdat_sp);
# dev.off()
#####################################################
## End: Generate a bubble plot suitable for slides
#####################################################


