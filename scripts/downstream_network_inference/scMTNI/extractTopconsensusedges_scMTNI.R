#!/user/bin/env Rscript

## Adapted from Shilu's script:
## /mnt/dv/wid/projects5/Roy-singlecell/shilu_work/gammlab/scripts/extractTopconsensusedges.R
##

args=commandArgs(trailingOnly = T);
if(length(args)==0)
{
  stop("At least one argument must be supplied(inputfile).n",call. = FALSE)
} else if(length(args)>0)
{
  print("Start!")
}

library(data.table)

# cl=args[1] ## cl="C1"
#n=args[2] ## n=cf0.5 or top5000
path=args[1] 
cf=args[2] #"top" or "cf"
n=as.numeric(args[3])
cellfile=args[4]

if (length(args) > 4) {
  datadir <- args[5] # contains ${cell}_genes_filteredlowexpression.txt. 
  # Note that for filtering out _TFA regulators with low expression, we will remove _TFA genes in ${cell}_genes_filteredlowexpression.txt (as they are activity inferred)
  # and then duplicate the remaining genes and appending _TFA to the duplicates. We remove _TFA regulators not within the list.
}

cells=read.table(cellfile,stringsAsFactors = F)
cells=cells$V1
data=NULL
for(cell in cells){
  # d=fread(paste0(path,cell,"/consensus_edges.txt"))
  d=fread(paste0(path,cell,"/consensus_edges_filteredlowexpression.txt"))
  d <- d[!grepl("_TFA$", V2)] # remove edges whose targets are _TFA features

  if (length(args) > 4) {
      geneskp <- fread(paste0(datadir, "/", cell, "_genes_filteredlowexpression.txt"), header=FALSE)
      geneskp <- geneskp$V1
      geneskp <- geneskp[!grepl("_TFA$", geneskp)] # Remove _TFA genes (as they are activity) from the filtered low expression list
      geneskp_TFA <- paste0(geneskp, "_TFA")
      d <- d[V1 %in% c(geneskp, geneskp_TFA), ]
  }

  d=d[order(d$V3,decreasing = T),]
  if(cf=="top"){

    ## Select the top n rows in d if n >= Number of rows in d;
    ## Otherwise, select all the rows in d.
    # d1=d[1:n,]
    ## d1=d[1:min(n, nrow(d)),] # only top 1000 not all edges with the same value, increases alphabetical bias - Suvojit
     top_n = min(n, nrow(d)) # Get the top n rows first
     threshold = d$V3[top_n] # value of the last edge in top n
     d1 = d[d$V3 >= threshold, ] # Select all edges >= that threshold
    ## alphabetical bias on top# edges is removed now - Suvojit dated October 31, 2025 during Brandlab analysis
  }else{ ## cf == "cf"
    id=which(d$V3>=n)
    d1=d[id,]
  }
  print(paste(cell,nrow(d1)))
  # fwrite(d1,file=paste0(path,cell,"/consensus_edges_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = F)
  fwrite(d1,file=paste0(path,cell,"/consensus_edges_filteredlowexpression_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = F)

  #d1=fread(paste0(path,cell,"/consensus_edges_top",n,".txt"))
  # print(paste0(path,cell,"/consensus_edges_",cf,n,".txt"))
  print(paste0(path,cell,"/consensus_edges_filteredlowexpression_",cf,n,".txt"))

  ## Set the col headers
  names(d1)=c("regulator","target",cell)

  if(is.null(data)){
    data=d1
  }else{
    data=merge(data,d1,by=c("regulator","target"),all=T)
  }
}
data[is.na(data)]=0
data$edge=paste0(data$regulator,"->",data$target)
selcols=c("edge",cells)
dall=data[,..selcols]
# fwrite(dall,file=paste0(path,"consensus_edges_merged_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = T)
fwrite(dall,file=paste0(path,"consensus_edges_merged_filteredlowexpression_",cf,n,".txt"),quote=F,sep="\t",row.names = F,col.names = T)


