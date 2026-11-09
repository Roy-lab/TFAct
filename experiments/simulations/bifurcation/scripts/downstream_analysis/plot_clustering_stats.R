library(ggplot2)
 
clustering_stat_file <- "TFA_output/analysis/clustering_stats.txt"

F_tab <-read.table(clustering_stat_file, header=TRUE, sep = "\t", check.names=FALSE)

F_tab$algo <- gsub("^.*/([^/]*)$","\\1", F_tab$algo)

F_tab$algo <- sapply(F_tab$algo, function(x){if(x=="hals_prior"){"HALS-p"}else{x}})

# ignore nm4f
F_tab <- F_tab[!F_tab$algo %in% c('ibpg', 'ibpg_wogr', 'ibpg_noreg', 'nmqf_hals', 'nmqf_hals_noreg'),]

system('mkdir -p TFA_output/analysis/clustering_stats')

ggplot(F_tab, aes(fill=algo, y=ARI, x=interaction(edge_addition_prob, lda_topics))) + 
    geom_boxplot(position=position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), text = element_text(size = 18))+
  coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(expand = c(0, -0.8))
ggsave("TFA_output/analysis/clustering_stats/ARI.png",width = 15, height = 10)


