library(ggplot2)
library(parallel)
library(tidyr)
library("scmamp")

#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/GRN_inference_w_TFA

theme_set(
    theme_minimal(base_family = "Arial")
)


# GRN_outputs/algos/inferelator_default/processed_FBS_reprogramming/cluster1/stability/mESC_chipunion_KDUnion_intersect/stats_10000/AUPR.txt
indir="GRN_outputs/algos"
fscores_tops = c(500, 1000, 2000)
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
    dat[[column_name]][dat[[column_name]] == "regularized_nca_tfa_and_exp"] <- "regularized_nca"
    dat[[column_name]][dat[[column_name]] == "regularized_nca_tfa_only"] <- "r_nca_tfa_only"
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
    c[c == "regularized_nca_tfa_only"] <- "r_nca_tfa_only"
    return(c)
}

datasets = c('processed_FBS_reprogramming', 'processed_Buenrostro_Hematopoiesis', 'processed_Ranzoni_fetal_hematopoiesis_fine_tree')
# algos=c('kl_nmf_2fct','HALS-p', 'ibpg', 'NINCA_quad', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2fct', 'ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'regularized_nca')
# algos=c('kl_nmf_2fct','kl_nmf_2fct_2000_1000','HALS-p', 'ibpg', 'NINCA_quad', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2fct', 'ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'regularized_nca')
algos=c('default', 'NINCA_quad', 'TIGER', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2fct', 'ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'no_TFA', 'kl_nmf_2fct','kl_nmf_2fct_2000_1000', 'NetREX_CF', 'regularized_nca_tfa_only', 'regularized_nca_tfa_and_exp')

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

system(paste0('mkdir -p ',indir,'/analysis/inferelator_AUPR'))
system(paste0('mkdir -p ',indir,'/analysis/inferelator_Fscore'))

all_mats <- list()
all_mats_Fscores <- list()
for(t in fscores_tops){
    all_mats_Fscores[[as.character(t)]] = list()
}

for(d in datasets){
    clusters <- cluster_lst[[d]]
    gold_standards <- gs_lst[[d]]
    cur_ds_mats <- list()
    cur_ds_Fscores_mats <- list()
    for(t in fscores_tops){
        cur_ds_Fscores_mats[[as.character(t)]] = list()
    }
    for(g in gold_standards){
        cur_mat <- matrix(0, length(clusters), length(algos))
        cur_mat_fs <- list()
        for(t in fscores_tops){
            cur_mat_fs[[as.character(t)]] = matrix(0, length(clusters), length(algos))
        }
        for(c in 1:length(clusters)){
            for(a in 1:length(algos)){
                cur_f <- paste0(indir,"/inferelator_",algos[[a]], '/',d,'/',clusters[[c]],'/stability/', g,'/stats_10000/AUPR.txt')
                cur_mat[c, a] = scan(cur_f, quiet = TRUE)
                for(t in fscores_tops){
                    cur_mat_fs[[as.character(t)]][c,a] = scan(paste0(indir,"/inferelator_",algos[[a]], '/',d,'/',clusters[[c]],'/stability/', g,'/stats_', t,'/stat.txt'), what = numeric(), sep = "\t", skip = 1, nlines = 1, quiet = TRUE)[1]
                }
            }
        }
        cur_mat <- as.data.frame(cur_mat)
        rownames(cur_mat) <- clusters
        colnames(cur_mat) <- algos
        if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree"){
            cur_mat$cell = rownames(cur_mat)
        } else {
            cur_mat$cell = gsub("^cluster", "C", rownames(cur_mat))
        }
        df <- pivot_longer(cur_mat, cols = -cell, names_to = "algorithm", values_to = "AUPR")
        df$AUPR <- as.numeric(df$AUPR)
        df <- formal_algo_names(df, "algorithm")
        midpoint <- (min(df$AUPR, na.rm = TRUE) + max(df$AUPR, na.rm = TRUE)) / 2
        df$text_color <- ifelse(df$AUPR > midpoint, "white", "black")
        ggplot(df, aes(x = cell, y = algorithm, fill = AUPR)) +
            geom_tile(color = "white") +  # White borders around tiles
            scale_fill_gradient(low = "white", high = "blue", limits = range(df$AUPR, na.rm = TRUE)) +  # Set max to 1
            geom_text(aes(label = round(AUPR, 3), color = text_color), size = 3) +  # Add numbers
            scale_color_identity() +
            labs(title = "", x = "", y = "", fill = "AUPR") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))  # Tilt x-axis labels
        if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree") {
            ggsave(paste0(indir,'/analysis/inferelator_AUPR/',d,'_',g,'.pdf'), width = 8, height = length(unique(df$algorithm)) * 0.3)
        } else {
            ggsave(paste0(indir,'/analysis/inferelator_AUPR/',d,'_',g,'.pdf'), width = 7, height = length(unique(df$algorithm)) * 0.3)
        }
        for(t in fscores_tops){
            cur_mat_fs[[as.character(t)]] <- as.data.frame(cur_mat_fs[[as.character(t)]])
            if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree"){
                rownames(cur_mat_fs[[as.character(t)]]) <- clusters
            } else {
                rownames(cur_mat_fs[[as.character(t)]]) = gsub("^cluster", "C", clusters)
            }
            colnames(cur_mat_fs[[as.character(t)]]) <- algos
            cur_mat_fs[[as.character(t)]]$cell = rownames(cur_mat_fs[[as.character(t)]])
            df <- pivot_longer(cur_mat_fs[[as.character(t)]], cols = -cell, names_to = "algorithm", values_to = "Fscore")
            df$Fscore <- as.numeric(df$Fscore)
            df <- formal_algo_names(df, "algorithm")
            midpoint <- (min(df$Fscore, na.rm = TRUE) + max(df$Fscore, na.rm = TRUE)) / 2
            df$text_color <- ifelse(df$Fscore > midpoint, "white", "black")
            ggplot(df, aes(x = cell, y = algorithm, fill = Fscore)) +
                geom_tile(color = "white") +  # White borders around tiles
                scale_fill_gradient(low = "white", high = "blue", limits = range(df$Fscore, na.rm = TRUE)) +  # Set max to 1
                geom_text(aes(label = round(Fscore, 3), color = text_color), size = 3) +  # Add numbers
                scale_color_identity() +
                labs(title = "", x = "", y = "", fill = "Fscore") +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 15))  # Tilt x-axis labels
            if(d == 'processed_Ranzoni_fetal_hematopoiesis_fine_tree'){
                ggsave(paste0(indir,'/analysis/inferelator_Fscore/',d,'_',g,'_', t, '.pdf'), width = 8, height = length(unique(df$algorithm)) * 0.3)
            } else {
                ggsave(paste0(indir,'/analysis/inferelator_Fscore/',d,'_',g,'_', t, '.pdf'), width = 7, height = length(unique(df$algorithm)) * 0.3)
            }

            cur_mat_fs[[as.character(t)]][, "cell"] <- NULL
            colnames(cur_mat_fs[[as.character(t)]]) <- col_formal_algo_names(colnames(cur_mat_fs[[as.character(t)]]))
            cur_ds_Fscores_mats[[as.character(t)]][[length(cur_ds_Fscores_mats[[as.character(t)]])+1]] <- cur_mat_fs[[as.character(t)]]

            print(paste0(g, " ", d, " ","Fscore top", "_", t, "Friedman Test: ", friedmanTest(cur_mat_fs[[as.character(t)]])))
            test <- nemenyiTest(cur_mat_fs[[as.character(t)]], alpha=0.05)
            # test$diff.matrix
            # abs(test$diff.matrix) > test$statistic

            # mats_df[,"regularized_nca_tfa_only"] <- NULL # unfair comparison as it has fewer tfs
            pdf(paste0(indir,'/analysis/inferelator_Fscore/Fscore_top_',t,'_',g,'_',d,'_critical_difference.pdf'), width = 10, height = 6)
            plotCD(cur_mat_fs[[as.character(t)]], alpha=0.05, cex=0.5, lwd = 2)
            dev.off()
        }
        # Friedman. https://rdrr.io/cran/scmamp/f/vignettes/Statistical_assessment_of_the_differences.Rmd
        cur_mat[,"cell"] <- NULL
        colnames(cur_mat) <- col_formal_algo_names(colnames(cur_mat))
        cur_ds_mats[[length(cur_ds_mats)+1]] <- cur_mat

        print(paste0(g, " ", d, " ", "AUPR Friedman Test: ", friedmanTest(cur_mat)))
        test <- nemenyiTest (cur_mat, alpha=0.05)
        # test$diff.matrix
        # abs(test$diff.matrix) > test$statistic
        pdf(paste0(indir,'/analysis/inferelator_AUPR/AUPR_',g,'_',d,'_critical_difference.pdf'), width = 10, height = 6)
        plotCD(cur_mat, alpha=0.05, cex=0.5, lwd = 2)
        dev.off()
    }
    cur_ds_mats_df <- do.call(rbind, cur_ds_mats)
    all_mats[[length(all_mats)+1]] <- cur_ds_mats_df
    # friedmanTest(cur_ds_mats_df)
    # test <- nemenyiTest(cur_ds_mats_df, alpha=0.05)
    # test$diff.matrix
    # abs(test$diff.matrix) > test$statistic
    # png('test.png', res=300, width = 3000, height = 1500)
    # plotCD(cur_ds_mats_df, alpha=0.05, cex=0.5, lwd = 2)
    # dev.off()
    for(t in fscores_tops){
        cur_ds_mats_Fscores_df <- do.call(rbind, cur_ds_Fscores_mats[[as.character(t)]])
        all_mats_Fscores[[as.character(t)]][[length(all_mats_Fscores[[as.character(t)]])+1]] = cur_ds_mats_Fscores_df
    }
}


all_mats_df <- do.call(rbind, all_mats)
print(paste0("AUPR overall Friedman Test: ", friedmanTest(all_mats_df)))
test <- nemenyiTest(all_mats_df, alpha=0.05)
# test$diff.matrix
# abs(test$diff.matrix) > test$statistic

# all_mats_df[,"regularized_nca_tfa_only"] <- NULL # unfair comparison as it has fewer tfs
pdf(paste0(indir,'/analysis/inferelator_AUPR/AUPR_critical_difference.pdf'), width = 10, height = 6)
plotCD(all_mats_df, alpha=0.05, cex=0.5, lwd = 2)
dev.off()

for (t in fscores_tops){
    mats_df <- do.call(rbind, all_mats_Fscores[[as.character(t)]])
    print(paste0("Fscore top", "_", t, "Friedman Test: ", friedmanTest(mats_df)))
    test <- nemenyiTest(mats_df, alpha=0.05)
    # test$diff.matrix
    # abs(test$diff.matrix) > test$statistic

    # mats_df[,"regularized_nca_tfa_only"] <- NULL # unfair comparison as it has fewer tfs
    pdf(paste0(indir,'/analysis/inferelator_Fscore/Fscore_top_',t,'_critical_difference.pdf'), width = 10, height = 6)
    plotCD(mats_df, alpha=0.05, cex=0.5, lwd = 2)
    dev.off()
}