## Input arguments:
## prefix = Title of the output network figure.
## data1 = Network edgeset information.
## outname = Name of the output network figure file.
## w = Desired width of the network figure file.
## h = Desired height of the network figure file.
##
## nr = Number of rows in the network figure. Each figure is expected to
## have multiple panels where each panel shows the giant components of a
## specific cell clusters. "nr" helps to decide how many panels will be
## there in each row.
##
## ntop = The number of regulators to be labeled in each cluster-specific
## network figure. If there are one or more regulators that have the
## same number of targets as the top "ntop"-th regulator, all these
## regulators will be labeled. Other regulators will not be labeled in
## the network figure.
##
## cellorder = The order of the cell clusters in which the figure panels
## should be ordered.
##
genNetall_multi=function(prefix,data1,outname,w,h,nr=2,ntop=20){
  # https://journal.r-project.org/archive/2017/RJ-2017-023/RJ-2017-023.pdf
  # https://cran.r-project.org/web/packages/ggnetwork/vignettes/ggnetwork.html#geom_edges
  library(dplyr)
  library(GGally)
  library(network)
  library(sna)
  library(ggplot2)
  library(tidyverse)
  library(gridExtra)
  library(ggnetwork)

  names(data1)=c("source","destination","weight","cell") ## edges
  data1$source=as.character(data1$source)
  data1$destination=as.character(data1$destination)
  # remove self loops
  data1 <- data1 %>% filter(source!=destination)
  ## get nodes:
  nodes=unique(c(data1$source,data1$destination))
  dnodes=data.frame(id=1:length(nodes),label=nodes) #nodesize=1
  edges <- data1 %>% left_join(dnodes, by = c("source" = "label")) %>% rename(from = id)
  edges <- edges %>% left_join(dnodes, by = c("destination" = "label")) %>% rename(to = id)
  edges <- select(edges, from, to, weight,cell)
  print(nrow(edges))
  dedges=unique(edges)
  dnodes$label=gsub("\\.","-",dnodes$label)
  mynetwork <- network(dedges, vertex.attr = dnodes, matrix.type = "edgelist", ignore.eval = FALSE,directed = T,multiple=T)
  set.edge.attribute(mynetwork, "edgePresense", ifelse(mynetwork %e% "weight" > 0, "edge present", "edge absent"))
  set.seed(10052016)
  gnet=ggnetwork(mynetwork, arrow.size = 0.1,arrow.gap = 0.015, by = "cell",weights = "weight",
                 layout = "fruchtermanreingold")

  ## get out-degree for each regulator node per cell:
  cells=unique(data1$cell)
  dnode=NULL
  for(cell in cells){
    id=which(data1$cell==cell)
    dcell=data1[id,]

    ## "dreg" is a data frame with two columns -- "source" and "ntarget."
    ## It lists all regulators in the giant components of the given cell cluster.
    ## The regulators are ordered in the descending order of their number of targets.
    dreg=dcell%>%group_by(source) %>% summarise(ntarget=n()) %>% arrange(desc(ntarget)) %>% as.data.frame()

    ## Assign labels to the regulators in "dreg."
    ## These labels will be printed on the network figure of the given cell cluster.
    ## Only top "ntop" number of regulators in "dreg" are assigned labels.
    ## If there are one or more regulators that have the same number of targets as
    ## the top "ntop"-th regulator, then these regulators are also assigned labels.
    ## The other regulators are not assigned labels. As a result, they will not be
    ## labeled in the network figure.
    dreg$labelshort=as.character(dreg$source)
    if(nrow(dreg)>ntop){
      cf=dreg$ntarget[ntop]
      dreg$labelshort[which(dreg$ntarget<cf)]=""
    }else if(nrow(dreg)==0){
      next
    }
    dreg$cell=cell
    dnode=rbind(dnode,dreg)
  }
  names(dnode)[1]="label"
  gnet=merge(gnet,dnode,by=c("label","cell"),all.x=T)
  gnet$ntarget[which(is.na(gnet$ntarget))]=0
  gnet$labelshort[which(is.na(gnet$labelshort))]=""  
  gnet$cell=gsub("cluster","C",gnet$cell)
  gnet$cell=factor(gnet$cell,levels = gsub("cluster","C",cells))
  cols <- c("edge present" = "salmon", "edge absent" = "grey")

  ## plot networks:
  #pdf(outname,width = 7,height = 6,useDingbats=FALSE)#pdf("test_MEF.pdf",width = 15,height = 10)
  png(outname,width = w,height = h,res=300)
  p=ggplot(gnet, aes(x, y, xend = xend, yend = yend)) +
    geom_edges(
      aes(color = edgePresense,alpha=weight),size=0.5, #size=weight+0.1
      arrow = arrow(length = unit(3, "pt"), type = "closed")) +
    geom_nodes(aes(size=ntarget),color = "lightskyblue",alpha=0.9) +  #aes(color = "lightskyblue"), 
    geom_nodetext(aes(label = labelshort,size=ntarget-1),color = "black")+  #fontface = "bold"
    theme(aspect.ratio=1,legend.position = "bottom")+
    scale_color_manual(values=cols)+ggtitle(prefix)+theme_void()+
  facet_wrap(~cell, nrow = nr) 
  #theme_facet(legend.position = "bottom",strip.text.x = element_text(size = 30,face="bold"))
  plot(p)
  dev.off()
}

