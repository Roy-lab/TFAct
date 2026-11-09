## TFA version of plot_topic_networks.R for LDA
## Adapted from Shilu's script:
## /mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/scripts/plotNetworks_LDA.R
##

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
genNetall_multi_path <- args[1]

cvnpath <- args[2] # "/mnt/dv/wid/projects7/Roy-singlecell2/suvo_work/scMTNI_brandlab_human_Combined/inputs_for_GRN_without_TFA/results_b_4_q_1/analysis/"

cellfile <- args[3] # /mnt/dv/wid/projects7/Roy-singlecell2/spencer_work/Brand_TAL1/brand_combine_order.txt

filter=args[4] # _cf0.7, _top5000

k <- as.numeric(args[5]) # number of topics

outdir <- args[6] # paste0(cvnpath,"lda_TFcellbygene/network",filter,"/k",k,"/")

if (length(args) > 6) {
  datadir <- args[7] # contains ${cell}_genes_filteredlowexpression.txt. 
  # Note that for filtering out _TFA regulators with low expression, we will remove _TFA genes in ${cell}_genes_filteredlowexpression.txt (as they are activity inferred)
  # and then duplicate the remaining genes and appending _TFA to the duplicates. We remove _TFA regulators not within the list.
  network_type <- args[8] # if "TFA", then filter for edges with _TFA regulators in the list above. If "gene", then filter for edges with gene regulators in the list above.
}

system(paste0("mkdir -p ", outdir))

indir = paste0(cvnpath,"/lda_TFcellbygene/network",filter,"/k",k,"/")

# /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/estimate_tfa/brand_multiome_project/scMTNI_with_TFA/scripts/genNetall_multi.R
source(genNetall_multi_path)

# cvnpath="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/Results/scCVN_upto1transition/liger_sqrt_ncell50_k8_filterhumanbc_FBS_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/"
# cvnpath="/mnt/dv/wid/projects5/Roy-singlecell/saptarshi_work/scMTNImanuscript/results/scMTNI/a2s/condorOutputs/analysis/"

# cellfile="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/liger_sqrt_ncell50_k10_filterhumanbc/celltype_order.txt"


# cells=read.table("/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/liger_sqrt_ncell50_k8_filterhumanbc_FBS/celltype_order.txt",stringsAsFactors = F)
cells=read.table(cellfile,stringsAsFactors = F)

cells=cells$V1

# filter="_full"



# ## get most divergent topic:
# dstat=NULL
# for(tp in 1:k){
#   for(cell in cells){
#     tfs=read.table(paste0(indir,"/TFs_topicid_",cell,".txt"))
#     tf=tfs[which(tfs$V2==tp),]
#     prefix=paste0("lda_topic",tp,"_",cell)
#     genes=read.table(paste0(indir,"/genes_topicid.txt"))
#     gene=genes[which(genes$V2==tp),]
#     dstat=rbind(dstat,data.frame(tp,cell,tf=nrow(tf),gene=nrow(gene)))
#   }
# }
# dt=dstat %>% group_by(tp) %>% summarise(mean=mean(tf),variance=var(tf)) %>% arrange(desc(variance))
# write.csv(dt,file=paste0(indir,"/lda_k",k,"_variance.csv"),quote=F,row.names=F)

## plot networks for each topic:
genes=read.table(paste0(indir,"/genes_topicid.txt"))
names(genes)=c("target","targetid","targetprob")

d_by_cell <- list()

for (cell in cells) {
  d <- fread(paste0(cvnpath, "/", cell, "/consensus_edges_filteredlowexpression.txt"))	
  names(d) <- c("regulator", "target", "weight")
  ## Remove TFA targets
  d <- d[!grepl("_TFA$", target), ]
  if (length(args) > 6) {
      geneskp <- fread(paste0(datadir, "/", cell, "_genes_filteredlowexpression.txt"), header=FALSE)
      geneskp <- geneskp$V1
      geneskp <- geneskp[!grepl("_TFA$", geneskp)] # Remove _TFA genes from the filtered low expression list
      geneskp_TFA <- paste0(geneskp, "_TFA")
    if (network_type == "TFA") { # TFA only edges
      d <- d[regulator %in% geneskp_TFA, ]
    } else if (network_type == "gene") { # gene only edges
      d <- d[regulator %in% geneskp, ]
    }
  }
  ## Apply filter
  if (grepl("^_top[0-9]+$", filter)) {
    top_n <- as.numeric(
      sub("^_top([0-9]+)$", "\\1", filter, perl = TRUE)
    )
    d <- d[order(-weight)]
    cutoff <- d$weight[top_n]
    d <- d[weight >= cutoff]
  } else if (grepl("^_cf0\\.[0-9]+$", filter)) {
    cutoff <- as.numeric(
      sub("^_cf([.0-9]+)$", "\\1", filter, perl = TRUE)
    )
    d <- d[weight >= cutoff]
  } else {
    stop("Error: filter must be either _top[0-9]+ or _cf0.[0-9]+")
  }
  d_by_cell[[cell]] <- d
}


for(tp in 1:k){ #c(2,3,10)
  ## extract networks for regulators and targets for each topic:
  dtopic=NULL
  gene=genes[which(genes$targetid==tp),]
  for(cell in cells){
    tfs=read.table(paste0(indir,"/TFs_topicid_",cell,".txt"))
    names(tfs)=c("regulator","tfid","tfprob")
    tf=tfs[which(tfs$tfid==tp),]
    print(paste(cell,"TF:",unique(tf$tfid),"gene:",unique(gene$targetid)))
    prefix=paste0("lda_topic",tp,"_",cell)

    d = d_by_cell[[cell]] # fread(paste0(cvnpath,cell,"/consensus_edges_filteredlowexpression",filter,".txt"))

    names(d)=c("regulator","target","weight")
    print(nrow(d))
    d2=merge(d,tfs,by="regulator")
    d1=merge(d2,genes,by="target")
        

    #id=which(d1$regulator%in%tf$regulator & d1$target%in%gene$target)
    id=which(d1$tfid==tp & d1$targetid==tp)
    d2=d1[id,]
    d2$cell=cell
    dtopic=rbind(dtopic,d2)
    
  }
  dtopic1=dtopic[,c("regulator","target","weight","cell")]
  fwrite(dtopic1,file=paste0(outdir,"/networks_lda_k",k,"_topic",tp,".txt"),quote=F,row.names=F,col.names=T,sep="\t")
  prefix=paste0("lda_topic",tp,filter)
  outname=paste0(outdir, '/',prefix,"_consensus_edges_filteredlowexpression_k",k,"_network.png")

  # options(error = recover)
  genNetall_multi(prefix,dtopic1,outname,w=5000,h=3000,nr=2,ntop=20)
}

