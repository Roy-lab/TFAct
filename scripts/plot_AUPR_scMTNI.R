library(ggplot2)
library(parallel)
library(tidyr)

#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data

# TIGER/{{dataset}}/{{c}}_stability/{{gs}}/stats_{{top_e}}/AUPR.txt", outdir=config["outdir"]),
indir="TFA_output"

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

datasets = c('processed_FBS_reprogramming', 'processed_Buenrostro_Hematopoiesis', 'processed_Ranzoni_fetal_hematopoiesis_fine_tree')
### for comparing against HALS-p and IBPG, use the definition of algos here
# algos=c('kl_nmf_2fct','kl_nmf_2fct_2000_1000','HALS-p', 'ibpg', 'NINCA_quad', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2fct', 'ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'regularized_nca')
#### for seeing the influce of hyperparameter choice on some other methods' performance, use this:
### also need to modify the ggsave below
# algos = c('NetREX_CF', 'NetREX_CF_1', 'NetREX_CF_90', 'TIGER_2_2_1_1_10', 'TIGER_1_1_2_2_10', 'TIGER_1_1_1_1_5')
#### otherwise:
algos=c('NINCA_quad', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2fct', 'ibpg_2fct_2000_1000', 'ibpg_2fct_2000_100', 'ibpg_2fct_100_1000', 'regularized_nca', 'prior')
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

system(paste0('mkdir -p ',indir,'/analysis/AUPR'))
system(paste0('mkdir -p ',indir,'/analysis/AUPR_Other_algos_hyperparam'))

for(d in datasets){
    clusters <- cluster_lst[[d]]
    gold_standards <- gs_lst[[d]]
    for(g in gold_standards){
        cur_mat <- matrix(0, length(clusters), length(algos))
        for(c in 1:length(clusters)){
            for(a in 1:length(algos)){
                cur_f <- paste0(indir,"/",algos[[a]], '/',d,'/',clusters[[c]],'_stability/', g,'/stats_10000/AUPR.txt')
                cur_mat[c, a] = scan(cur_f, quiet = TRUE)
            }
        }
        cur_mat <- as.data.frame(cur_mat)
        if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree"){
            rownames(cur_mat) <- clusters
        } else {
            rownames(cur_mat) <- gsub("^cluster", "C", clusters)
        }
        colnames(cur_mat) <- algos
        cur_mat$cell = rownames(cur_mat)
        df <- pivot_longer(cur_mat, cols = -cell, names_to = "algorithm", values_to = "AUPR")
        df$AUPR <- as.numeric(df$AUPR)
        df <- formal_algo_names(df, "algorithm")
        # df$text_color <- ifelse(df$AUPR > 0.5, "white", "black")
        midpoint <- (min(df$AUPR, na.rm = TRUE) + max(df$AUPR, na.rm = TRUE)) / 2
        df$text_color <- ifelse(df$AUPR > midpoint, "white", "black")
        ggplot(df, aes(x = cell, y = algorithm, fill = AUPR)) +
            geom_tile(color = "white") +  # White borders around tiles
            scale_fill_gradient(low = "white", high = "blue", limits = range(df$AUPR, na.rm = TRUE)) +  # Set max to 1
            geom_text(aes(label = round(AUPR, 2), color = text_color),  size = 6) +  # Add numbers
            scale_color_identity() +
            labs(title = "", x = "", y = "", fill = "AUPR") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 28))  # Tilt x-axis labels
        # ggsave(paste0(indir,'/analysis/AUPR/',d,'_',g,'.png'))
        if (d == "processed_Ranzoni_fetal_hematopoiesis_fine_tree") {
            ggsave(paste0(indir,'/analysis/AUPR/',d,'_',g,'.pdf'), width = 12, height = 7)
        } else {
            ggsave(paste0(indir,'/analysis/AUPR/',d,'_',g,'.pdf'), width = 10, height = 6.5)
        }
        # ggsave(paste0(indir,'/analysis/AUPR_Other_algos_hyperparam/',d,'_',g,'.png')) # use this for TIGER and NetREX_CF's impact of hyperparameter choice on performance
    }
}