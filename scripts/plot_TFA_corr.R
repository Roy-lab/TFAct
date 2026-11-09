library(ggplot2)
library(parallel)

# TODO: for all runs
formal_algo_names <- function(dat, column_name){
    dat[[column_name]][dat[[column_name]] == "ibpg"] <- "IBPG-p"
    dat[[column_name]][dat[[column_name]] == "hals_prior"] <- "HALS-p"
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

algos=c('baselines/hals_prior', 'ibpg', 'ibpg_wogr', 'ibpg_noreg', 'nmqf_hals', 'nmqf_hals_noreg', 'FastNCA', 'NINCA_quad', 'pmf_grn', 'regularized_nca', 'BITFAM', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2f', 'ibpg_2fct', 'kl_nmf_2fct')
ignored_algos <- c('ibpg_2f', 'ibpg_noreg', 'ibpg_wogr', 'nmqf_hals', 'nmqf_hals_noreg', 'ibpg', 'hals_prior')
# ignored_algos <- NULL
#algos=c('baselines/hals_prior', 'FastNCA', 'NINCA_quad', 'pmf_grn', 'regularized_nca', 'BITFAM', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2f')
Edge_adding_probs=c(0.2, 0.4)
Ldaks=c(10,20)
Prefixes=c('Spom', 'Anc2', 'Anc5', 'Calb', 'Cgla', 'right', 'Sbay')
# Exp_mat_w_dropouts=c('E', 'E_0.3', 'E_0.6')
Exp_mat_w_dropouts=c('E', 'E_0.6')
TFA_sim_reps=c('rep0')

result<- expand.grid(Edge_adding_probs, Ldaks, Prefixes, Exp_mat_w_dropouts, TFA_sim_reps)
list_of_vectors <- split(result, seq(nrow(result)))

plot_TFA <- function(edge_adding_prob, ldak, prefix, exp_mat, sim_rep, outdir){

    combined_df <-data.frame()

    for(i in 1:length(algos)){
        cur_f <- paste0("TFA_output/",algos[i],"/",edge_adding_prob,"/k",ldak,"/",prefix,"/rep0/",exp_mat,"_stability/stats_1000/tf_specifics.txt") # note that stats_1000, stats_2000's tf_specifics are the same as the number just affects fscores.
        cur_df <-read.table(cur_f, header=TRUE, sep = "\t", check.names=FALSE)

        cur_df$algo <- gsub("^.*/([^/]*)$","\\1",algos[i])

        combined_df <- rbind(combined_df, cur_df)
    }
    combined_df <- combined_df[!combined_df$algo %in% ignored_algos, ]
    combined_df <- formal_algo_names(combined_df, "algo")

    # suffix=paste0(edge_adding_prob, "_k", ldak,"_",prefix,"_",sim_rep,"_",exp_mat,"_1000.png")
    suffix=paste0(edge_adding_prob, "_k", ldak,"_",prefix,"_",sim_rep,"_",exp_mat,"_1000.pdf")

    ggplot(combined_df, aes(y=Pearson, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25), text = element_text(size = 25))

    ggsave(paste0("TFA_output/analysis/tfa_corr/Pearson_",suffix))
    ggplot(combined_df, aes(y=Spearman, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25), text = element_text(size = 25))

    ggsave(paste0("TFA_output/analysis/tfa_corr/Spearman_",suffix))

    return(combined_df[,c('algo', 'Pearson', 'Spearman')])
}

plot_all_TFA_by_group <- function(tab, suffix){
    tab <- tab[!tab$algo %in% ignored_algos, ]
    tab <- formal_algo_names(tab, 'algo')
    ggplot(tab, aes(y=Pearson, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25), text = element_text(size = 25)) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_cartesian(ylim = c(-1, 1))

    ggsave(paste0("TFA_output/analysis/tfa_corr/Pearson_",suffix))

    ggplot(tab, aes(y=Spearman, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 25), text = element_text(size = 25)) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + coord_cartesian(ylim = c(-1, 1))

    ggsave(paste0("TFA_output/analysis/tfa_corr/Spearman_",suffix))
}

outdir='TFA_output/analysis/tfa_corr'
system(paste0('mkdir -p ', outdir))

stat_lst<-mclapply(list_of_vectors, function(x){plot_TFA(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]])}, mc.cores = 20)

# covariates that may affect TFA
cov_classes <- lapply(list_of_vectors, function(x){paste(x$Var1, x$Var2, x$Var4, sep=",")}) # Var1 is 0.2/0.4, Var2 is 20/10, Var4 is expression, potentially with dropouts

unique_classes <- unique(cov_classes)

mclapply(c(1:length(unique_classes)), function(x){plot_all_TFA_by_group(do.call(rbind, stat_lst[cov_classes %in% unique_classes[[x]]]), paste0(gsub(",", "_", unique_classes[[x]]), ".pdf"))})
# mclapply(c(1:length(unique_classes)), function(x){plot_all_TFA_by_group(do.call(rbind, stat_lst[cov_classes %in% unique_classes[[x]]]), paste0(gsub(",", "_", unique_classes[[x]]), ".png"))})