library(ggplot2)
library(parallel)
 
Fscore_file <- "TFA_output/analysis/ct/Fscores.txt"

F_tab <-read.table(Fscore_file, header=TRUE, sep = "\t", check.names=FALSE)

F_tab$algo <- gsub("^.*/([^/]*)$","\\1", F_tab$algo)

F_tab$algo <- paste(F_tab$algo, F_tab$l_prior, F_tab$l_sim, sep="_")

formal_algo_names <- function(dat, column_name){
  for(i in c(0, 10, 100, 1000)){
    for(j in c(0, 10, 100, 1000)){
      dat[[column_name]][dat[[column_name]] == paste("ibpg_2fct", i, j,sep="_")] <- paste("TFAct-F", i, j,sep="_")
      dat[[column_name]][dat[[column_name]] == paste("kl_nmf_2fct", i, j,sep="_")] <- paste("TFAct-KL", i, j,sep="_")
    }
  }   
  return(dat)
}

F_tab <- formal_algo_names(F_tab, "algo")

F_tab <- F_tab[!F_tab$exp_dropout=="E_0.3",]
F_tab$exp_dropout <- gsub('E$', 'nodrop',F_tab$exp_dropout)
F_tab$exp_dropout <- gsub('E_', 'dropout',F_tab$exp_dropout)
# for(E in unique(F_tab$exp_dropout)){
#   cur_F_tab <- F_tab[F_tab$exp_dropout == E,]
#   ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_gt, x=interaction(Fscore_top_edge, lda_topics))) +
#       geom_boxplot(position=position_dodge(width = 0.8)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 18, face="bold"))+
#     coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(expand = c(0, -0.8))
#   ggsave(paste0("TFA_output/analysis/ct/Fscores_ground_true", gsub("^E", "", E), ".png"), width = 15, height = 10)

#   ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_prior, x=interaction(Fscore_top_edge, lda_topics))) +
#       geom_boxplot(position=position_dodge(width = 0.8)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 18, face="bold"))+
#     coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(expand = c(0, -0.8))
#   ggsave(paste0("TFA_output/analysis/ct/Fscores_prior", gsub("^E", "", E), ".png"), width = 15, height = 10)
# }



dummy <- mclapply(unique(F_tab$Fscore_top_edge), function(E){
                for(ad in unique(F_tab$edge_addition_prob)){
                  cur_F_tab <- F_tab[F_tab$Fscore_top_edge == E & F_tab$edge_addition_prob == ad,]
                  ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_gt, x=interaction(exp_dropout, lda_topics))) +
                      geom_boxplot(position=position_dodge(width = 0.8), outlier.shape = NA) +
                    theme(text = element_text(size = 50))+ 
                    scale_x_discrete(expand = c(0, -0.8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_cartesian(ylim = c(0, 0.51))
                  # ggsave(paste0("TFA_output/analysis/ct/Fscores_ground_true",  E, "_", ad, ".png"), width = 15, height = 10) 
                  ggsave(paste0("TFA_output/analysis/ct/Fscores_ground_true",  E, "_", ad, ".pdf"), width = 25, height = 10) 


                  ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_prior, x=interaction(exp_dropout, lda_topics))) +
                      geom_boxplot(position=position_dodge(width = 0.8), outlier.shape = NA) + 
                    theme(text = element_text(size = 50))+
                    scale_x_discrete(expand = c(0, -0.8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_cartesian(ylim = c(0, 0.51))
                  # ggsave(paste0("TFA_output/analysis/ct/Fscores_prior", E, "_", ad, ".png"), width = 15, height = 10)
                  ggsave(paste0("TFA_output/analysis/ct/Fscores_prior", E, "_", ad, ".pdf"), width = 25, height = 10)
                }
    }, mc.cores = 50)