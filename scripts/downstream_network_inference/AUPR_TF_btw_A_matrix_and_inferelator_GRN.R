
A_mat_indir = '../scMTNI_data/TFA_output/' # resolved
Inferelator_indir = 'GRN_outputs/algos'
system('mkdir -p GRN_outputs/algos/analysis/inferelator_AUPR_TF/comparison_w_Amat')
system('mkdir -p GRN_outputs/algos/analysis/inferelator_AUPR_TF/comparison_w_Amat_overleaf_friendly')

datasets = c('processed_FBS_reprogramming', 'processed_Buenrostro_Hematopoiesis', 'processed_Ranzoni_fetal_hematopoiesis_fine_tree')

mouse_gold_standards = c('mESC_chipunion_KDUnion_intersect', 'mESC_KDUnion', 'mESC_chipunion')

human_gold_standards = c('Cus_KO_0.01_Cus_ChIP_2_intersect.gs_network', 'Cus_ChIP_2.gs_network', 'Cus_KO_0.01.gs_network', 'Cus_KO_0.01_union_intersect.gs_network', 'UniBind.gs_network')

gs_lst = list(mouse_gold_standards, human_gold_standards, human_gold_standards)
names(gs_lst) <- datasets

all_A_mats <- list()
all_inferelator_mats <- list()
pvals_lst <- list()

for(d in datasets){
    for(g in gs_lst[[d]]){
        for(stat in c('median', 'mean')){
            cur_a_mat <- read.table(paste0(A_mat_indir,'/analysis/AUPR_TF/', stat, '_',d,'_',g,'_tab.txt'), sep = "\t", header = TRUE,
                                row.names = 1, check.names = FALSE)
            cur_inferelator <- read.table(paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/', stat, '_',d,'_',g,'_tab.txt'), sep = "\t", header = TRUE,
                                row.names = 1, check.names = FALSE)
            cur_inferelator$regularized_nca <- NULL
            inferelator_colnames <- colnames(cur_inferelator)
            inferelator_colnames[inferelator_colnames == "regularized_nca"] <- "regularized_nca_w_expression" # for convention
            inferelator_colnames[inferelator_colnames == "r_nca_tfa_only"] <- "regularized_nca" # as A matrix did the filtering for regularized nca
            colnames(cur_inferelator) <- inferelator_colnames
            intersected_method = intersect(colnames(cur_a_mat), colnames(cur_inferelator)) 
            cur_a_mat = cur_a_mat[, intersected_method]
            cur_inferelator = cur_inferelator[, intersected_method]

            all_A_mats[[d]][[g]][[stat]] = cur_a_mat
            all_inferelator_mats[[d]][[g]][[stat]] = cur_inferelator

            # test the hypothesis that performance increases
            pvals <- sapply(intersected_method, function(m){wilcox.test(cur_a_mat[, m], cur_inferelator[, m], paired = TRUE, alternative = "less")$p.value})
            pvals_df <- data.frame(algorithm = names(pvals), pvalue = as.numeric(pvals))
            write.table(pvals_df, paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/comparison_w_Amat/', stat, '_',d,'_',g,'_wilcox_signed_rank.txt'), sep = "\t",
                row.names = FALSE, quote = FALSE)
            
            pvals_lst[[d]][[g]] = pvals_df
            colnames(pvals_lst[[d]][[g]]) <- c('Algorithm', gsub('.gs_network$', '', g))
            rownames(pvals_lst[[d]][[g]]) <- pvals_lst[[d]][[g]][['Algorithm']]
            pvals_lst[[d]][[g]][['Algorithm']] <- NULL
            # pvals_df$pvalue = paste0(pvals_df$pvalue, '\\\\')
            # pvals_df$algorithm <- gsub('_', '\\\\_', pvals_df$algorithm)
            # write.table(pvals_df, paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/comparison_w_Amat_overleaf_friendly/', stat, '_',d,'_',g,'_wilcox_signed_rank.txt'), sep = " & ",
            #     row.names = FALSE, quote = FALSE)
        }
    }
}

lapply(datasets, function(d){
    cur_df <- Reduce(function(x, y){
        y = y[rownames(x), , drop = FALSE]
        return(cbind(x,y))
    }, pvals_lst[[d]])
    colnames(cur_df) <- gsub('_', '\\\\_', colnames(cur_df))
    cur_df[[length(cur_df)]] <- paste0(cur_df[[length(cur_df)]], '\\\\')
    rownames(cur_df) <- gsub('_', '\\\\_', rownames(cur_df))
    write.table(cur_df, paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/comparison_w_Amat_overleaf_friendly/combined_', stat, '_',d, '_wilcox_signed_rank.txt'), sep = " & ",
            row.names = TRUE, quote = FALSE, col.names=NA)
})


# per dataset
for(d in datasets){
    for(stat in c('median', 'mean')){
        cur_ds_gs = gs_lst[[d]]
        cur_a_mat <- do.call(rbind, lapply(cur_ds_gs, function(g){all_A_mats[[d]][[g]][[stat]]}))
        cur_inferelator <- do.call(rbind, lapply(cur_ds_gs, function(g){all_inferelator_mats[[d]][[g]][[stat]]}))

        pvals <- sapply(colnames(cur_a_mat), function(m){wilcox.test(cur_a_mat[, m], cur_inferelator[, m], paired = TRUE, alternative = "less")$p.value})
        pvals_df <- data.frame(algorithm = names(pvals), pvalue = as.numeric(pvals))
        write.table(pvals_df, paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/comparison_w_Amat/', stat, '_',d,'_wilcox_signed_rank.txt'), sep = "\t",
            row.names = FALSE, quote = FALSE)
        
        pvals_df$pvalue = paste0(pvals_df$pvalue, '\\\\')
        pvals_df$algorithm <- gsub('_', '\\\\_', pvals_df$algorithm)
        write.table(pvals_df, paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/comparison_w_Amat_overleaf_friendly/', stat, '_',d,'_all_gold_standards_wilcox_signed_rank.txt'), sep = " & ",
                row.names = FALSE, quote = FALSE)
    }
}
# per gold standard (human only)
for(g in human_gold_standards){
    for(stat in c('median', 'mean')){
        cur_g_ds <- c('processed_Buenrostro_Hematopoiesis', 'processed_Ranzoni_fetal_hematopoiesis_fine_tree')
        cur_a_mat <- do.call(rbind, lapply(cur_g_ds, function(d){all_A_mats[[d]][[g]][[stat]]}))
        cur_inferelator <- do.call(rbind, lapply(cur_g_ds, function(d){all_inferelator_mats[[d]][[g]][[stat]]}))

        pvals <- sapply(colnames(cur_a_mat), function(m){wilcox.test(cur_a_mat[, m], cur_inferelator[, m], paired = TRUE, alternative = "less")$p.value})
        pvals_df <- data.frame(algorithm = names(pvals), pvalue = as.numeric(pvals))
        write.table(pvals_df, paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/comparison_w_Amat/', stat, '_',g,'_wilcox_signed_rank.txt'), sep = "\t",
            row.names = FALSE, quote = FALSE)
        # pvals_df$pvalue = paste0(pvals_df$pvalue, '\\\\')
        # pvals_df$algorithm <- gsub('_', '\\\\_', pvals_df$algorithm)
        # write.table(pvals_df, paste0(Inferelator_indir,'/analysis/inferelator_AUPR_TF/comparison_w_Amat_overleaf_friendly/', stat, '_',g,'_wilcox_signed_rank.txt'), sep = " & ",
        #         row.names = FALSE, quote = FALSE)
    }
}




            