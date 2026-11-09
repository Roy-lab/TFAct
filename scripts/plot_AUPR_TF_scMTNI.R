library(ggplot2)
library(ggrastr)
library(parallel)
library(tidyr)
library("scmamp")

theme_set(
    theme_minimal(base_family = "Arial")
)

#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data

# TIGER/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
indir="TFA_output"
max_cores=10

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

col_formal_algo_names <- function(c){
    c[c == "ibpg_2fct"] <- "TFAct-F_500_100"
    c[c == "ibpg_2fct_500_1000"] <- "TFAct-F_500_1000"
    c[c == "ibpg_2fct_500_0"] <- "TFAct-F_500_0"
    c[c == "ibpg_2fct_2000_100"] <- "TFAct-F_2000_100"
    c[c == "ibpg_2fct_2000_1000"] <- "TFAct-F_2000_1000"
    c[c == "ibpg_2fct_100_1000"] <- "TFAct-F_100_1000"
    c[c == "kl_nmf_2fct"] <- "TFAct-KL_500_100"
    c[c == "kl_nmf_2fct_2000_1000"] <- "TFAct-KL_2000_1000"
    c[c == "regularized_nca_tfa_and_exp"] <- "regularized_nca"
    return(c)
}

datasets = c('processed_FBS_reprogramming', 'processed_Buenrostro_Hematopoiesis', 'processed_Ranzoni_fetal_hematopoiesis_fine_tree')
# algos=c('kl_nmf_2fct','kl_nmf_2fct_2000_1000','HALS-p', 'ibpg', 'NINCA_quad', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2fct', 'ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'regularized_nca')
algos=c('NINCA_quad', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2fct', 'ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'regularized_nca', 'prior', "kl_nmf_2fct", "kl_nmf_2fct_2000_1000")
AUPR_TF_dir = 'AUPR_TF'

# #### for seeing the influce of hyperparameter choice on some other methods' performance, use this:
# algos = c('NetREX_CF', 'NetREX_CF_1', 'NetREX_CF_90', 'TIGER_2_2_1_1_10', 'TIGER_1_1_2_2_10', 'TIGER_1_1_1_1_5')
# AUPR_TF_dir = 'AUPR_TF_Other_algos_hyperparam'

mouse_gold_standards = c('mESC_chipunion_KDUnion_intersect', 'mESC_KDUnion', 'mESC_chipunion')

human_gold_standards = c('Cus_KO_0.01_Cus_ChIP_2_intersect.gs_network', 'Cus_ChIP_2.gs_network', 'Cus_KO_0.01.gs_network', 'Cus_KO_0.01_union_intersect.gs_network', 'UniBind.gs_network')

# FBS reprogramming
FBS_clusters = c("cluster1","cluster2","cluster3","cluster4","cluster5","cluster7","cluster8")

# Buenrostro_Hematopoiesis 
BH_clusters = c("cluster1","cluster2","cluster3","cluster6","cluster7","cluster8", "cluster9", "cluster10")

# Ranzoni_fetal_hematopoiesis/coarse_tree
RFH_fine_clusters = c("B-cells","Erythroid-cells","GPs","Granulocytes","HSC-MPPs","LMPs","Mast-cells","Megakaryocytes","MEMPs","Monocytes","pDCs")

cluster_lst = list(FBS_clusters, BH_clusters, RFH_fine_clusters)
names(cluster_lst) <- datasets
gs_lst = list(mouse_gold_standards, human_gold_standards, human_gold_standards)
names(gs_lst) <- datasets

system(paste0('mkdir -p ',indir,'/analysis/', AUPR_TF_dir))

all_mats_median <- list()
all_mats_mean <- list()
all_mats_per_TF <- list()
for(d in datasets){
    clusters <- cluster_lst[[d]]
    gold_standards <- gs_lst[[d]]
    cur_ds_mats <- mclapply(gold_standards, function(g){
        cur_mat <- matrix(0, length(clusters), length(algos))
        cur_mat_mean <- matrix(0, length(clusters), length(algos))
        all_tf_focused_mat_cur_c <- list()
        for(c in 1:length(clusters)){
            tf_focused_mat <- NULL
            for(a in 1:length(algos)){
                cur_f <- paste0(indir,"/",algos[[a]], '/',d,'/',clusters[[c]],'_stability/', g,'/stats_10000/AUPR_TF.txt')
                cur_f <- read.table(cur_f,sep='\t', header=FALSE) 
                cur_mat[c, a] = median(cur_f$V2)
                cur_mat_mean[c, a] = mean(cur_f$V2)
                colnames(cur_f) <- c('TF', algos[[a]])
                if (algos[[a]] == "regularized_nca"){ # this method trims regulators.
                    next
                }
                if (is.null(tf_focused_mat)){
                    tf_focused_mat <- cur_f
                } else {
                    tf_focused_mat <- merge(tf_focused_mat, cur_f, by='TF', all=FALSE)
                }
            }
            all_tf_focused_mat_cur_c[[length(all_tf_focused_mat_cur_c)+1]] = tf_focused_mat # for each cluster, we treat each TF as a comparison
        }
        if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree"){
            names(all_tf_focused_mat_cur_c) <- clusters
        } else {
            names(all_tf_focused_mat_cur_c) <- gsub("^cluster", "C", clusters)
        }
        all_tf_focused_mat_cur_c_df <- do.call(rbind, Map(function(df, name) {
            df$cluster <- name
            df
        }, all_tf_focused_mat_cur_c, names(all_tf_focused_mat_cur_c)))

        long_all_tf_focused_mat_df <- pivot_longer(all_tf_focused_mat_cur_c_df[, setdiff(names(all_tf_focused_mat_cur_c_df), c('TF'))],
                   cols = -cluster,
                    names_to = "algo",
                    values_to = "AUPR_TF"
                    )
        long_all_tf_focused_mat_df <- formal_algo_names(long_all_tf_focused_mat_df, "algo")
        ggplot(long_all_tf_focused_mat_df, aes(x = algo, y = AUPR_TF, color = cluster)) +
            geom_point_rast(
                position = position_jitterdodge(
                        jitter.width = 0.1,
                        dodge.width = 0.8
                    ),
                size = 1.5,alpha = 0.6, raster.dpi = 300) +
            labs(x = "", y = "", title = "") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 25),panel.spacing = unit(0, "lines"))  
        ggsave(paste0(indir,'/analysis/', AUPR_TF_dir, '/per_TF_',d,'_',g,'.pdf'), height = 8, width = length(unique(long_all_tf_focused_mat_df$algo)) * 0.7)

        all_tf_focused_mat_cur_c_df[, "cluster"] <- NULL

        cur_mat <- as.data.frame(cur_mat)
        if( d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree"){
            rownames(cur_mat) <- clusters
        } else {
            rownames(cur_mat) <- gsub("^cluster", "C", clusters)
        }
        colnames(cur_mat) <- algos
        cur_mat$cell = rownames(cur_mat)
        df <- pivot_longer(cur_mat, cols = -cell, names_to = "algorithm", values_to = "AUPR")
        df$AUPR <- as.numeric(df$AUPR)
        df <- formal_algo_names(df, "algorithm")
        midpoint <- (min(df$AUPR, na.rm = TRUE) + max(df$AUPR, na.rm = TRUE)) / 2
        df$text_color <- ifelse(df$AUPR > midpoint, "white", "black")
        ggplot(df, aes(x = cell, y = algorithm, fill = AUPR)) +
            geom_tile(color = "white") +  # White borders around tiles
            scale_fill_gradient(low = "white", high = "blue", limits = range(df$AUPR, na.rm = TRUE)) +  # Set max to 1
            geom_text(aes(label = round(AUPR, 2), color = text_color), size = 3) +  # Add numbers
            scale_color_identity() +
            labs(title = "", x = "", y = "", fill = "medianAUPR") +
            theme_minimal()  +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15),
          panel.spacing = unit(0, "lines"))  
        if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree") {
            ggsave(paste0(indir,'/analysis/', AUPR_TF_dir, '/median_',d,'_',g,'.pdf'), width = 7, height = length(unique(df$algorithm)) * 0.3)
        } else {
            ggsave(paste0(indir,'/analysis/', AUPR_TF_dir, '/median_',d,'_',g,'.pdf'), width = 6, height = length(unique(df$algorithm)) * 0.3)
        }
        cur_mat_mean <- as.data.frame(cur_mat_mean)
        if( d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree"){
            rownames(cur_mat_mean) <- clusters
        } else {
            rownames(cur_mat_mean) <- gsub("^cluster", "C", clusters)
        }
        colnames(cur_mat_mean) <- algos
        cur_mat_mean$cell = rownames(cur_mat_mean)
        df <- pivot_longer(cur_mat_mean, cols = -cell, names_to = "algorithm", values_to = "AUPR")
        df$AUPR <- as.numeric(df$AUPR)
        df <- formal_algo_names(df, "algorithm")
        midpoint <- (min(df$AUPR, na.rm = TRUE) + max(df$AUPR, na.rm = TRUE)) / 2
        df$text_color <- ifelse(df$AUPR > midpoint, "white", "black")
        ggplot(df, aes(x = cell, y = algorithm, fill = AUPR)) +
            geom_tile(color = "white") +  # White borders around tiles
            scale_fill_gradient(low = "white", high = "blue", limits = range(df$AUPR, na.rm = TRUE)) +  # Set max to 1
            geom_text(aes(label = round(AUPR, 2), color = text_color), size = 3) +  # Add numbers
            scale_color_identity() +
            labs(title = "", x = "", y = "", fill = "meanAUPR") +
            theme_minimal()  +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15),
          panel.spacing = unit(0, "lines"))  
        if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree") {
            ggsave(paste0(indir,'/analysis/', AUPR_TF_dir, '/mean_',d,'_',g,'.pdf'), width = 7, height = length(unique(df$algorithm)) * 0.3)
        } else {
            ggsave(paste0(indir,'/analysis/', AUPR_TF_dir, '/mean_',d,'_',g,'.pdf'), width = 6, height = length(unique(df$algorithm)) * 0.3)
        }
        
        cur_mat[, "cell"] <- NULL
        cur_mat_mean[, "cell"] <- NULL
        colnames(cur_mat) <- col_formal_algo_names(colnames(cur_mat))
        colnames(cur_mat_mean) <- col_formal_algo_names(colnames(cur_mat_mean))

        all_tf_focused_mat_cur_c_df_2 <- all_tf_focused_mat_cur_c_df
        all_tf_focused_mat_cur_c_df_2[,"TF"] <- NULL
        colnames(all_tf_focused_mat_cur_c_df_2) <- col_formal_algo_names(colnames(all_tf_focused_mat_cur_c_df_2))

        # friedmanTest(all_mats_df)
        # test <- nemenyiTest(all_mats_df, alpha=0.05)
        # all_mats_df[,"regularized_nca_tfa_only"] <- NULL # unfair comparison as it has fewer tfs
        pdf(paste0(indir,'/analysis/', AUPR_TF_dir, '/median_',d,'_',g,'_critical_difference.pdf'), width = 10, height = 6)
        plotCD(cur_mat, alpha=0.05, cex=0.5, lwd = 2)
        dev.off()

        pdf(paste0(indir,'/analysis/', AUPR_TF_dir, '/mean_',d,'_',g,'_critical_difference.pdf'), width = 10, height = 6)
        plotCD(cur_mat_mean, alpha=0.05, cex=0.5, lwd = 2)
        dev.off()

        pdf(paste0(indir,'/analysis/',AUPR_TF_dir,'/per_TF_',d,'_',g,'_critical_difference.pdf'), width = 10, height = 6)
        plotCD(all_tf_focused_mat_cur_c_df_2, alpha=0.05, cex=0.5, lwd = 2)
        dev.off()

        write.table(cur_mat, paste0(indir,'/analysis/',AUPR_TF_dir,'/median_',d,'_',g,'_tab.txt'), sep = "\t", quote = FALSE, col.names = NA)
        write.table(cur_mat_mean, paste0(indir,'/analysis/',AUPR_TF_dir,'/mean_',d,'_',g,'_tab.txt'), sep = "\t", quote = FALSE, col.names = NA)

        return(list(cur_mat, cur_mat_mean, all_tf_focused_mat_cur_c_df))
    }, mc.cores=max_cores)
    all_mats_median[[length(all_mats_median)+1]] = do.call(rbind, lapply(cur_ds_mats, function(cur_i){cur_i[[1]]}))
    all_mats_mean[[length(all_mats_mean)+1]] = do.call(rbind, lapply(cur_ds_mats, function(cur_i){cur_i[[2]]}))
    all_mats_per_TF[[length(all_mats_per_TF)+1]] = do.call(rbind, lapply(cur_ds_mats, function(cur_i){cur_i[[3]]}))
}

# median
all_mats_df <- do.call(rbind, all_mats_median)
friedmanTest(all_mats_df)
test <- nemenyiTest(all_mats_df, alpha=0.05)
# test$diff.matrix
# abs(test$diff.matrix) > test$statistic

# all_mats_df[,"regularized_nca_tfa_only"] <- NULL # unfair comparison as it has fewer tfs
pdf(paste0(indir,'/analysis/',AUPR_TF_dir,'/median_critical_difference.pdf'), width = 10, height = 6)
plotCD(all_mats_df, alpha=0.05, cex=0.5, lwd = 2)
dev.off()

# mean
all_mats_df <- do.call(rbind, all_mats_mean)
friedmanTest(all_mats_df)
test <- nemenyiTest(all_mats_df, alpha=0.05)
# test$diff.matrix
# abs(test$diff.matrix) > test$statistic

# all_mats_df[,"regularized_nca_tfa_only"] <- NULL # unfair comparison as it has fewer tfs
pdf(paste0(indir,'/analysis/',AUPR_TF_dir,'/mean_critical_difference.pdf'), width = 10, height = 6)
plotCD(all_mats_df, alpha=0.05, cex=0.5, lwd = 2)
dev.off()

# per TF
all_mats_df <- do.call(rbind, all_mats_per_TF)
all_mats_df[,"TF"] <- NULL
colnames(all_mats_df) <- col_formal_algo_names(colnames(all_mats_df))
friedmanTest(all_mats_df)
test <- nemenyiTest(all_mats_df, alpha=0.05)
# test$diff.matrix
# abs(test$diff.matrix) > test$statistic

# all_mats_df[,"regularized_nca_tfa_only"] <- NULL # unfair comparison as it has fewer tfs
pdf(paste0(indir,'/analysis/',AUPR_TF_dir,'/per_TF_critical_difference.pdf'), width = 10, height = 6)
plotCD(all_mats_df, alpha=0.05, cex=0.5, lwd = 2)
dev.off()