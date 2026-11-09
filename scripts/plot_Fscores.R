library(ggplot2)
library(parallel)

# theme_set(
#     theme_minimal(base_family = "Arial")
# )
 
Fscore_file <- "TFA_output/analysis/Fscores.txt"

F_tab <-read.table(Fscore_file, header=TRUE, sep = "\t", check.names=FALSE)

F_tab$algo <- gsub("^.*/([^/]*)$","\\1", F_tab$algo)

F_tab$algo <- sapply(F_tab$algo, function(x){if(x=="hals_prior"){"HALS-p"}else{x}})

formal_algo_names <- function(dat, column_name){
    dat[[column_name]][dat[[column_name]] == "ibpg"] <- "IBPG-p"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct"] <- "TFAct-F_500_100"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_1000"] <- "TFAct-F_500_1000"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_0"] <- "TFAct-F_500_0"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_2000_100"] <- "TFAct-F_2000_100"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_2000_1000"] <- "TFAct-F_2000_1000"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_100_1000"] <- "TFAct-F_100_1000"
    dat[[column_name]][dat[[column_name]] == "kl_nmf_2fct"] <- "TFAct-KL_500_100"
    dat[[column_name]][dat[[column_name]] == "kl_nmf_2fct_2000_1000"] <- "TFAct-KL_2000_1000"
    return(dat)
}

# ignore nm4f
# F_tab <- F_tab[!F_tab$algo %in% c('ibpg', 'ibpg_wogr', 'ibpg_noreg', 'nmqf_hals', 'nmqf_hals_noreg'),]


system(paste0('mkdir -p ', 'TFA_output/analysis/gt'))
system(paste0('mkdir -p ', 'TFA_output/analysis/prior'))
###
## things commented out below were used to plot the old plot where every top edges were included in the same plot
###

# for(E in unique(F_tab$exp_dropout)){
#   cur_F_tab <- F_tab[F_tab$exp_dropout == E,]
#   ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_gt, x=interaction(Fscore_top_edge, lda_topics))) +
#       geom_boxplot(position=position_dodge(width = 0.8)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))+
#     coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(expand = c(0, -0.8))
#   ggsave(paste0("TFA_output/analysis/Fscores_ground_true", gsub("^E", "", E), ".png"), width = 15, height = 10)

#   ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_prior, x=interaction(Fscore_top_edge, lda_topics))) +
#       geom_boxplot(position=position_dodge(width = 0.8)) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))+
#     coord_cartesian(ylim = c(0, 1)) + scale_x_discrete(expand = c(0, -0.8))
#   ggsave(paste0("TFA_output/analysis/Fscores_prior", gsub("^E", "", E), ".png"), width = 15, height = 10)
# }

F_tab <- formal_algo_names(F_tab, "algo")
ignored_algos <- c('ibpg_2f', 'ibpg_noreg', 'ibpg_wogr', 'nmqf_hals', 'nmqf_hals_noreg', 'HALS-p', 'IBPG-p')
F_tab <- F_tab[!F_tab$algo %in% ignored_algos, ]
F_tab <- F_tab[!F_tab$exp_dropout=='E_0.3',]
F_tab$exp_dropout <- gsub('E$', 'nodrop',F_tab$exp_dropout)
F_tab$exp_dropout <- gsub('E_', 'dropout',F_tab$exp_dropout)

dummy <- mclapply(unique(F_tab$Fscore_top_edge), function(E){
                for(ad in unique(F_tab$edge_addition_prob)){
                  cur_F_tab <- F_tab[F_tab$Fscore_top_edge == E & F_tab$edge_addition_prob == ad,]
                  dropout_map <- c(
                      nodrop = "No Dropout",
                      dropout0.6 = "Dropout 0.6"
                    )
                  topic_map <- c(
                      k10 = "Higher noise",
                      k20 = "Lower noise"
                    )
                  if(ad == 0.2){
                    ymax = 0.75
                  } else {
                    ymax = 0.5
                  }
                  
                  ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_gt, x=interaction(dropout_map[as.character(exp_dropout)], topic_map[as.character(lda_topics)], sep = "\n"))) +
                      geom_boxplot(position=position_dodge(width = 0.8), outlier.shape = NA) +
                    theme(text = element_text(size = 30), axis.text.y = element_text(size = 40))+ 
                    scale_x_discrete(expand = c(0, -0.8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_cartesian(ylim = c(0, ymax))
                  # ggsave(paste0("TFA_output/analysis/gt/Fscores_top",  E, "_", ad, ".png"), width = 15, height = 10)
                  ggsave(paste0("TFA_output/analysis/gt/Fscores_top",  E, "_", ad, ".pdf"), width = 15, height = 10)

                  ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_prior, x=interaction(dropout_map[as.character(exp_dropout)], topic_map[as.character(lda_topics)], sep = "\n"))) +
                      geom_boxplot(position=position_dodge(width = 0.8), outlier.shape = NA) + 
                    theme(text = element_text(size = 30), axis.text.y = element_text(size = 40))+
                    scale_x_discrete(expand = c(0, -0.8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_cartesian(ylim = c(0, ymax))
                  # ggsave(paste0("TFA_output/analysis/prior/Fscores_top", E, "_", ad, ".png"), width = 15, height = 10)
                  ggsave(paste0("TFA_output/analysis/prior/Fscores_top", E, "_", ad, ".pdf"), width = 15, height = 10)
                }
    }, mc.cores = 50)