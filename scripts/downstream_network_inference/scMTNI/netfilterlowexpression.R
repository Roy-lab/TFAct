#!/user/bin/env Rscript
args=commandArgs(trailingOnly = T);
if(length(args)==0)
{
  stop("At least one argument must be supplied(inputfile).n",call. = FALSE)
}else if(length(args)>0)
{
  print("Start!")
}

cellfile=args[1]  ##/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/celltype_order.txt
cvnpath=args[2]  ##/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/Results/scCVN_upto1transition/liger_4cells_sqrt_19genesrm_varthre0.05_k10_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/
indir=args[3] ## directory for *.table ##"/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/"
if(length(args)==4){
  lowexprgene=1
}else{
  lowexprgene=0
}
library(data.table)
cells=read.table(cellfile,stringsAsFactors = F)
cells=cells$V1
#cvnpath="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/Results/scCVN_upto1transition/liger_4cells_sqrt_19genesrm_varthre0.05_k10_macs2/pg0.2_pm0.8_pr0.2_maxReg50_b4_bm4/subsample/analysis/"
#cells=read.table("/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/celltype_order.txt",stringsAsFactors = F)
filterlowexpression=function(indir,cell){
  #indir="/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/Buenrostro_2018/data/CVNdata/liger_4cells_sqrt_19genesrm_varthre0.05_k10/"
  dt=fread(paste0(indir,cell,".table")) 
  dt$Gene=gsub(paste0("_",cell),"",dt$Gene)
  rna=rowSums(dt[,-1]>0)
  df=data.frame(gene=dt$Gene,ncell=rna)
  idr=which(df$ncell<5)
  print(paste(cell,"removed",length(idr),"genes"))
  if(length(idr)>0){
  geneskeep=df$gene[-idr] } else {
  geneskeep=df$gene
  }
  write.table(geneskeep,file=paste0(indir,cell,"_genes_filteredlowexpression.txt"),quote=F,sep="\n",row.names = F,col.names =F)
  regulators=fread(paste0(indir,cell,"_allregulators.txt"),header=F)
  regulators$V1=gsub(paste0("_",cell),"",regulators$V1)
  idrg=which(regulators$V1%in%geneskeep)
  regulatorskeep=regulators$V1[idrg]
  write.table(regulatorskeep,file=paste0(indir,cell,"_regulators_filteredlowexpression.txt"),quote=F,sep="\n",row.names = F,col.names =F)
  return(data.frame(cell,geneskeep=length(geneskeep),regulatorskeep=length(regulatorskeep),ngene=nrow(dt),nreg=nrow(regulators)))
}

if(lowexprgene==1){
  ## filter low expression genes:
  dstat=NULL
  for(cell in cells){
    dd=filterlowexpression(indir,cell)
    dstat=rbind(dstat,dd)
  }
  print(dstat)
}

## get mean or median or 80% confident edges:
nedge=NULL
for(cell in cells){
  ## filter low expression genes:
  geneskp=fread(paste0(indir,cell,"_genes_filteredlowexpression.txt"),header=F)
  geneskp=geneskp$V1
  infile=paste0(cvnpath,cell,"/consensus_edges.txt")
  d=fread(infile)
  idk=which(d$V1%in%geneskp & d$V2%in%geneskp)
  print(paste("original net",nrow(d),"edges","removed",nrow(d)-length(idk),"edges"))
  d1=d[idk,]  
  d1=d1[order(d1$V3,decreasing=T),]
  fwrite(d1,file=paste0(cvnpath,cell,"/consensus_edges_filteredlowexpression.txt"),quote=F,sep="\t",row.names = F,col.names =F)
  print(paste0(cvnpath,cell,"/consensus_edges_filteredlowexpression.txt"))
  id=which(d1$V3>=0.8)
  nedge=c(nedge,nrow(d1[id,]))
}
print("80% confident edges")
print(paste(cells,nedge))
n=as.integer(max(nedge)/1000)*1000  
print(paste0("top edges: ",n))


