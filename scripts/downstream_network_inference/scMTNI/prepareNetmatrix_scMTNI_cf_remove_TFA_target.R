## Adapted from Sunnie Grace's script:
## /mnt/dv/wid/projects7/Roy-singlecell2/ashtonlab/results/H72_de/scmtni/merlin/scripts/prepareNetmatrix.R
##

#!/user/bin/env Rscript
args=trimws(commandArgs(trailingOnly = T));
if(length(args)==0)
{
  stop("At least one argument must be supplied(inputfile).n",call. = FALSE)
}else if(length(args)>0)
{
  print("Start!")
}

cellfile=args[1]  ##/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/celltype_order.txt
indir=args[2]  ##/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/Results/scCVN_upto1transition/liger_4cells_sqrt_19genesrm_varthre0.05_k10_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/
filter=args[3] #"_full" or "_cf0.8" or _top5000
if(length(args)==3){
  binary=""
}else{
  binary=args[4]
}

# print("Done reading in files")

library(data.table)
library(Matrix)
#cells=read.table("/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/celltype_order.txt",stringsAsFactors = F)
cells=read.table(cellfile,stringsAsFactors = F)
cells=cells$V1
#indir="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/Results/scCVN_upto1transition/liger_4cells_sqrt_19genesrm_varthre0.05_k10_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/"
#filter="_full"
#binary=""
#filter="_cf0.8"
#binary="_binary"
outdir=paste0(indir,"lda_TFcellbygene/network",filter,binary,"/")
dir.create(outdir,recursive = T)
print(outdir)
genes=NULL
regulators=NULL
for(cell in cells){
  #d=fread(paste0(indir,cell,"/consensus_edges_filteredlowexpression.txt")) #consensus_edges.txt
  # d=fread(paste0(indir,cell,"/consensus_edges",filter,".txt"))
  if ( grepl("^_top[0-9]+$", filter)){
    d=fread(paste0(indir,cell,"/consensus_edges_filteredlowexpression.txt"))	
    d = d[!grepl("_TFA$", d[[2]]), ]

    top_n = as.numeric(sub("^_top([.0-9]+)$", "\\1", filter, perl=TRUE))
    d = d[order(-d[[3]])]
    cutoff = d[[3]][top_n]
    d = d[d[[3]] >= cutoff]

  }
  else { #cf
    d=fread(paste0(indir,cell,"/consensus_edges_filteredlowexpression.txt"))	
    d = d[!grepl("_TFA$", d[[2]]), ]
    if(grepl("^_cf0.[0-9]+$", filter)){
          id=which(d$V3>=as.numeric(sub("^_cf([.0-9]+)$", "\\1", filter,perl=TRUE)))
        d=d[id,]
    } else {
      stop("Error: filter must be either _top[0-9]+ or _cf0.[0-9]+")
    }
  }
  regulators=c(regulators,unique(d$V1))
  reg_count=length(unique(d$V1))
  print(paste0("The number of regulators for ",cell," are ",reg_count))
  genes=c(genes,unique(d$V2))
  gene_count=length(unique(d$V2))
  print(paste0("The number of targets for ",cell," are ",gene_count))
}
genes=unique(genes)

genes=data.frame(gene=genes,id=1:length(genes))
genes$gene=as.character(genes$gene)
write.table(genes,file=paste0(outdir,"genes_filteredlowexpression",filter,binary,".txt"),quote=F,row.names=F,col.names=F,sep="\t")

regulators=sort(unique(regulators))
regulators=data.frame(gene=regulators, id = 1:length(regulators))
# regulators=merge(regulators,genes,by="gene")
write.table(regulators,file=paste0(outdir,"TFtarget_regulator_id",filter,binary,".txt"),quote=F,row.names=F,col.names=F,sep="\t")
# regulators$id=1:nrow(regulators)
regulators$gene=as.character(regulators$gene)
for(cell in cells){
  if ( grepl("^_top[0-9]+$", filter)){
    d=fread(paste0(indir,cell,"/consensus_edges_filteredlowexpression.txt"))	
    d = d[!grepl("_TFA$", d[[2]]), ]

    top_n = as.numeric(sub("^_top([.0-9]+)$", "\\1", filter, perl=TRUE))
    d = d[order(-d[[3]])]
    cutoff = d[[3]][top_n]
    d = d[d[[3]] >= cutoff]
  } else { #cf
    d=fread(paste0(indir,cell,"/consensus_edges_filteredlowexpression.txt"))
    d = d[!grepl("_TFA$", d[[2]]), ]
    if(grepl("^_cf0.[0-9]+$", filter)){
          id=which(d$V3>=as.numeric(sub("^_cf([.0-9]+)$", "\\1", filter,perl=TRUE)))
        d=d[id,]
    }
  }
  if(binary=="_binary"){
    d$V3=1
  }
  regs=data.frame(gene=unique(d$V1))
  regs=merge(regulators,regs,by="gene")
  regs=regs[order(regs$id),]
  #print(paste0("The number of regulators for ",cell," are ",regs,"\n"))
  #regs$id=1:nrow(regs)
  #dmat=matrix(0,nrow(regs),nrow(genes))
  dmat=matrix(0,nrow(regulators),nrow(genes))
  for(i in regs$id) #1:nrow(regs))
  {
    #reg=regs$gene[i]
    reg=regulators$gene[i]
    d1=d[which(d$V1==reg),]
    d1=merge(d1,genes,by.x="V2",by.y="gene")
    dmat[i,d1$id]=d1$V3
  }
  #dmat=forceSymmetric(dmat, "UP")
  dmat=data.table(as.matrix(dmat))
  row.names(dmat)=regulators$gene  #regs$gene 
  names(dmat)=genes$gene
  fwrite(dmat,file=paste0(outdir,cell,"_consensus_edges_filteredlowexpression_mat",filter,binary,".txt"),quote=F,row.names=T,col.names=T,sep="\t")
}


