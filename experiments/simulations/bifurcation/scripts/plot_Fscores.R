library(ggplot2)

 
Fscore_file <- "TFA_output/analysis/Fscores.txt"

F_tab <-read.table(Fscore_file, header=TRUE, sep = "\t", check.names=FALSE)

F_tab$algo <- gsub("^.*/([^/]*)$","\\1", F_tab$algo)

F_tab$algo <- sapply(F_tab$algo, function(x){if(x=="hals_prior"){"HALS-p"}else{x}})

# ignore nm4f
# F_tab <- F_tab[!F_tab$algo %in% c('ibpg', 'ibpg_wogr', 'ibpg_noreg', 'nmqf_hals', 'nmqf_hals_noreg'),]

ggplot(F_tab, aes(fill=algo, y=fscore_w_gt, x=interaction(Fscore_top_edge, lda_topics))) + 
    geom_boxplot(position=position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))+
  coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(expand = c(0, -0.8))
ggsave("TFA_output/analysis/Fscores_ground_true.pdf",width = 15, height = 10)
# ggsave("TFA_output/analysis/Fscores_ground_true.png",width = 15, height = 10)

ggplot(F_tab, aes(fill=algo, y=fscore_w_prior, x=interaction(Fscore_top_edge, lda_topics))) + 
    geom_boxplot(position=position_dodge(width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))+
  coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(expand = c(0, -0.8))
ggsave("TFA_output/analysis/Fscores_prior.pdf",width = 15, height = 10)
# ggsave("TFA_output/analysis/Fscores_prior.png",width = 15, height = 10)

