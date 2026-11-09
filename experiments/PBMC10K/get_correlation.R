library(dplyr)
library(parallel)
library(ggplot2)
library(tidyr)
max_cores=50

#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/PBMC10k
filtered_chrom_assay <- read.table('chromVAR_TFA.txt', header = TRUE,  sep = '\t',  stringsAsFactors = FALSE, check.names=FALSE)
all_reg <- read.table("TIGER_zenodo/data/PBMC_ct/all_regulators.txt")$V1
filtered_chrom_assay <- filtered_chrom_assay[all_reg, ] # reorder for matching P.txt


outdir="TFA_output"
celltypes=c("CD4_Naive","CD4_TCM","CD8_Naive","CD16_Mono","NK","Treg","CD14_Mono","cDC","CD8_TEM_1","Intermediate_B",
        "Memory_B","Plasma","CD4_TEM","MAIT","Naive_B","CD8_TEM_2","pDC","gdT","HSPC")
cellnames = lapply(celltypes, function(x){read.table(paste0('TIGER_zenodo/data/PBMC_ct/', x, '_cellnames.txt'))$V1})
names(cellnames) <- celltypes
stability_selection_algos=c('HALS-p', 'ibpg_2fct', 'ibpg', 'ibpg_2fct_500_1000', 'ibpg_2fct_500_0')
algo_w_TF_names=c('regularized_nca')
other_algos=c('TIGER','supirfactor_hierarchical', 'supirfactor_shallow', 'NetREX_CF','NINCA_quad') 
algos = c(stability_selection_algos, algo_w_TF_names, other_algos)

ct_gt_TFA <- lapply(celltypes, function(x){filtered_chrom_assay[, cellnames[[x]]]})
names(ct_gt_TFA) <- celltypes

getcorr <- function(c, alg){
    if (alg %in% stability_selection_algos){
        cur_P <- read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'_stability/P.txt'), check.names=FALSE, header=FALSE)
    } else if (alg %in% algo_w_TF_names){
        cur_P <- read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'_stability/P.txt'), check.names=FALSE, header=FALSE)
        rownames(cur_P) <- cur_P$V1
        cur_P$V1 <- NULL
        nonactive_tfs <- setdiff(all_reg, rownames(cur_P))
        nonactive_tfs <- as.data.frame(matrix(0, length(nonactive_tfs), ncol(cur_P)), row.names=nonactive_tfs)
        colnames(cur_P) <- colnames(nonactive_tfs) 
        cur_P <- rbind(cur_P, nonactive_tfs)
        cur_P <- cur_P[all_reg,]
    }
    else if(alg %in% other_algos){
        cur_P <- read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'/P.txt'), check.names=FALSE, header=FALSE)
    } else {
        stop("unsupport algorithm name")
    }
    # scale and threshold positive TFA values as some does not do things like l1 selection
    if (!(alg %in% c('NexREX_CF', 'FastNcA', 'NINCA_quad', 'regularized_nca', 'supirfactor_shallow'))){
        row_sums <- rowSums(abs(cur_P))
        cur_P[row_sums > 0, ] <- cur_P[row_sums > 0, ] / row_sums[row_sums > 0]
        cur_P[cur_P<0.5/ncol(cur_P)] = 0
    }
    cur_gt <- ct_gt_TFA[[c]]
    pearson_correlations <- sapply(seq(1:nrow(cur_P)), 
            function(i){cor(unlist(cur_P[i, ]), unlist(cur_gt[i, ]), method='pearson')}
        )
    Spearman_correlations <- sapply(seq(1:nrow(cur_P)), 
            function(i){cor(unlist(cur_P[i, ]), unlist(cur_gt[i, ]), method='spearman')}
        )
    return(list(pearson_correlations,Spearman_correlations))
}


corr_lst<-lapply(algos, function(a){mclapply(celltypes, function(c){getcorr(c, a)}, mc.cores=max_cores)})
names(corr_lst) <- algos

system('mkdir -p TFA_output/analysis/TFA_corr')
pearson_lst_ct <- list()
spearman_lst_ct <- list()

# grouped by algorithms, individual cell types
for(c in 1:length(celltypes)){
    pearson_lst_ct[[length(pearson_lst_ct)+1]] <- do.call(rbind, lapply(algos, function(a){corr_lst[[a]][[c]][[1]]}))
    spearman_lst_ct[[length(spearman_lst_ct)+1]] <- do.call(rbind,lapply(algos, function(a){corr_lst[[a]][[c]][[2]]}))
}

formal_algo_names <- function(dat, column_name){
    dat[[column_name]][dat[[column_name]] == "ibpg"] <- "IBPG-p"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct"] <- "TFAct-F_500_100"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_1000"] <- "TFAct-F_500_1000"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_0"] <- "TFAct-F_500_0"
    return(dat)
}

plot_c_corr <- function(c){
    cur_pearson <- as.data.frame(pearson_lst_ct[[c]])
    cur_pearson$algo <- algos
    cur_pearson <- pivot_longer(
                            cur_pearson,
                            cols = -algo,      
                            values_to = "pearson_c"         
                            )
    cur_pearson <- formal_algo_names(cur_pearson, "algo")
    ggplot(cur_pearson, aes(x = algo, y = pearson_c)) +
        geom_boxplot()  + ggtitle(paste0("Pearson correlation coefficients for ", celltypes[[1]])) + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
    ggsave(paste0('TFA_output/analysis/TFA_corr/Pearson_', celltypes[[c]],'.pdf'))

    ggplot(cur_pearson, aes(x = algo, y = pearson_c)) +
        stat_summary(fun = mean, geom = "bar", fill = "skyblue", color = "black", width = 0.7) +
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                    geom = "errorbar", width = 0.2, color = "red") +
        stat_summary(fun = mean, geom = "text", 
               aes(label = round(..y.., 4)), 
               vjust = -0.5, color = "black", size = 3.5) +
        labs(
            title = paste0("mean Pearson correlation coefficients for ", celltypes[[1]]),
            x = "algorithm",
            y = "Pearson cc"
        ) + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0('TFA_output/analysis/TFA_corr/mean_Pearson_', celltypes[[c]],'.pdf'))

    cur_spearman <- as.data.frame(spearman_lst_ct[[c]])
    cur_spearman$algo <- algos
    cur_spearman <- pivot_longer(
                            cur_spearman,
                            cols = -algo,      
                            values_to = "spearman_c"         
                            )
    
    cur_spearman <- formal_algo_names(cur_spearman, "algo")
    ggplot(cur_spearman, aes(x = algo, y = spearman_c)) +
        geom_boxplot()  + ggtitle(paste0("Spearman correlation coefficients for ", celltypes[[1]])) + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
    ggsave(paste0('TFA_output/analysis/TFA_corr/Spearman_', celltypes[[c]],'.pdf'))

        ggplot(cur_spearman, aes(x = algo, y = spearman_c)) +
        stat_summary(fun = mean, geom = "bar", fill = "skyblue", color = "black", width = 0.7) +
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                    geom = "errorbar", width = 0.2, color = "red") +
        stat_summary(fun = mean, geom = "text", 
               aes(label = round(..y.., 4)), 
               vjust = -0.5, color = "black", size = 5) +
        labs(
            title = paste0("mean Spearman correlation coefficients for ", celltypes[[1]]),
            x = "algorithm",
            y = "Spearman cc"
        ) + theme_minimal() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0('TFA_output/analysis/TFA_corr/mean_Spearman_', celltypes[[c]],'.pdf'))
}

mclapply(1:length(celltypes), function(c){plot_c_corr(c)}, mc.cores=max_cores)


# grouped by algorithms, collected over all cell types. This may be futile as # cells differ greatly among cell types
pearson_collapse <- as.data.frame(do.call(cbind, pearson_lst_ct))
pearson_collapse$algo <- algos
pearson_collapse <- pivot_longer(
                    pearson_collapse,
                    cols = -algo,      
                    values_to = "pearson_c"         
                    )
    ggplot(pearson_collapse, aes(x = algo, y = pearson_c)) +
        geom_boxplot()  + ggtitle("Pearson correlation coefficients, aggregated over cell types") + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
ggsave('TFA_output/analysis/TFA_corr/Pearson.pdf')

pearson_collapse_group_mean<- aggregate(x= pearson_collapse$pearson_c,
                      by = list(pearson_collapse$algo),      
                      FUN = function(x){mean(x, na.rm=TRUE)})
pearson_collapse_group_sd<- aggregate(x= pearson_collapse$pearson_c,
                      by = list(pearson_collapse$algo),      
                      FUN = function(x){sd(x, na.rm=TRUE)})

pearson_collapse_group <- data.frame(algo=pearson_collapse_group_mean[,1], mean=pearson_collapse_group_mean$x, sd=pearson_collapse_group_sd$x)


spearman_collapse <- as.data.frame(do.call(cbind, spearman_lst_ct))
spearman_collapse$algo <- algos
spearman_collapse <- pivot_longer(
                    spearman_collapse,
                    cols = -algo,      
                    values_to = "spearman_c"         
                    )
    ggplot(spearman_collapse, aes(x = algo, y = spearman_c)) +
        geom_boxplot()  + ggtitle("Spearman correlation coefficients, aggregated over cell types") + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
ggsave('TFA_output/analysis/TFA_corr/Spearman.pdf')

spearman_collapse_group_mean<- aggregate(x= spearman_collapse$spearman_c,
                      by = list(spearman_collapse$algo),      
                      FUN = function(x){mean(x, na.rm=TRUE)})
spearman_collapse_group_sd<- aggregate(x= spearman_collapse$spearman_c,
                      by = list(spearman_collapse$algo),      
                      FUN = function(x){sd(x, na.rm=TRUE)})

spearman_collapse_group <- data.frame(algo=spearman_collapse_group_mean[,1], mean=spearman_collapse_group_mean$x, sd=spearman_collapse_group_sd$x)


##
## Correlation across cell types. No scaling
##

get_concate_P <- function(alg) {
    if (alg %in% stability_selection_algos){
        cur_P <- do.call(cbind, mclapply(celltypes, function(c){read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'_stability/P.txt'), check.names=FALSE, header=FALSE)}, mc.cores=max_cores))
    } else if (alg %in% algo_w_TF_names){
        cur_P <- do.call(cbind, mclapply(celltypes, function(c){
            cur_P_ <- read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'_stability/P.txt'), check.names=FALSE, header=FALSE)
            rownames(cur_P_) <- cur_P_$V1
            cur_P_$V1 <- NULL
            nonactive_tfs <- setdiff(all_reg, rownames(cur_P_))
            nonactive_tfs <- as.data.frame(matrix(0, length(nonactive_tfs), ncol(cur_P_)), row.names=nonactive_tfs)
            colnames(cur_P_) <- colnames(nonactive_tfs) 
            cur_P_ <- rbind(cur_P_, nonactive_tfs)
            cur_P_ <- cur_P_[all_reg,]
            }, mc.cores=max_cores))
    }
    else if(alg %in% other_algos){
        cur_P <- do.call(cbind, mclapply(celltypes, function(c){read.table(paste0(outdir,'/',alg,'/PBMC_ct/',c,'/P.txt'), check.names=FALSE, header=FALSE)}, mc.cores=max_cores))
    } else {
        throw("unsupport algorithm name")
    }

    # does not scale and remove low TFA here, as TF with high scaled activity may be low in a global sense
    return(cur_P)
}

Ps_cellnames <- unlist(cellnames)
cell_reordered_filtered_chrom_assay <- filtered_chrom_assay[, Ps_cellnames]
system('mkdir -p TFA_output/analysis/TFA_corr/overall')


# takes a bit of time but doable in < 5 mins. 
pearson_lst <- rep(list(0), length(algos))
spearman_lst <- rep(list(0), length(algos))
names(pearson_lst) <- algos
names(spearman_lst) <- algos

results_corr <- mclapply(algos, function(a){
    cur_P <- get_concate_P(a)
    pearson_correlations <- sapply(seq(1:nrow(cur_P)), 
            function(i){cor(unlist(cur_P[i, ]), unlist(cell_reordered_filtered_chrom_assay[i, ]), method='pearson')}
        )
    write.table(pearson_correlations, paste0('TFA_output/analysis/TFA_corr/overall/Pearson_', a, '.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
    Spearman_correlations <- sapply(seq(1:nrow(cur_P)), 
            function(i){cor(unlist(cur_P[i, ]), unlist(cell_reordered_filtered_chrom_assay[i, ]), method='spearman')}
        )
    write.table(Spearman_correlations, paste0('TFA_output/analysis/TFA_corr/overall/Spearman_', a, '.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
    return(list(pearson = pearson_correlations, spearman = Spearman_correlations))
}, mc.cores=max_cores)

for (i in seq_along(algos)) {
    pearson_lst[[algos[i]]] <- results_corr[[i]]$pearson
    spearman_lst[[algos[i]]] <- results_corr[[i]]$spearman
}

pearson_tab <- unlist(pearson_lst)
pearson_tab <- data.frame(algo=unlist(lapply(algos, function(a){rep(a, length(all_reg))})), pearson=pearson_tab)

pearson_tab <-  formal_algo_names(pearson_tab, 'algo')

ggplot(pearson_tab, aes(x = algo, y = pearson)) +
        geom_boxplot()  + ggtitle("Pearson correlation coefficients over the full datasets") + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
ggsave('TFA_output/analysis/TFA_corr/overall/Pearson.pdf')

pearson_group_mean<- aggregate(x= pearson_tab$pearson,
                      by = list(pearson_tab$algo),      
                      FUN = function(x){mean(x, na.rm=TRUE)})
pearson_group_sd<- aggregate(x= pearson_tab$pearson,
                      by = list(pearson_tab$algo),      
                      FUN = function(x){sd(x, na.rm=TRUE)})
pearson_group <- data.frame(algo=pearson_group_mean[,1], mean=pearson_group_mean$x, sd=pearson_group_sd$x)


spearman_tab <- unlist(spearman_lst)
spearman_tab <- data.frame(algo=unlist(lapply(algos, function(a){rep(a, length(all_reg))})), spearman=spearman_tab)

spearman_tab <- formal_algo_names(spearman_tab, 'algo')

ggplot(spearman_tab, aes(x = algo, y = spearman)) +
        geom_boxplot()  + ggtitle("Spearman correlation coefficients over the full datasets") + theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1) )
ggsave('TFA_output/analysis/TFA_corr/overall/Spearman.pdf')

spearman_group_mean<- aggregate(x= spearman_tab$spearman,
                      by = list(spearman_tab$algo),      
                      FUN = function(x){mean(x, na.rm=TRUE)})
spearman_group_sd<- aggregate(x= spearman_tab$spearman,
                      by = list(spearman_tab$algo),      
                      FUN = function(x){sd(x, na.rm=TRUE)})
spearman_group <- data.frame(algo=spearman_group_mean[,1], mean=spearman_group_mean$x, sd=spearman_group_sd$x)

###
# 2. global mean plot
###

pearson_group <- formal_algo_names(pearson_group, "algo")
spearman_group <- formal_algo_names(spearman_group, "algo")

excluded_algo <- c('HALS-p', 'IBPG-p')

pearson_group <- pearson_group[!pearson_group$algo %in% excluded_algo, ]
spearman_group <- spearman_group[!spearman_group$algo %in% excluded_algo,]

ggplot(pearson_group, aes(x = algo, y = mean)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "red") +
  labs(
    title = "Global Pearson CC over the full datasets",
    x = "",
    y = ""
  ) + geom_text(aes(label = round(mean, 4)), vjust = -0.5, size = 5, color = "black") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14) )
ggsave("TFA_output/analysis/TFA_corr/overall/mean_Pearson.pdf")

ggplot(spearman_group, aes(x = algo, y = mean)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, color = "red") +
  labs(
    title = "Global Spearman CC over the full datasets",
    x = "",
    y = ""
  ) + geom_text(aes(label = round(mean, 4)), vjust = -0.5, size = 5, color = "black") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 14) )
ggsave("TFA_output/analysis/TFA_corr/overall/mean_Spearman.pdf")


###
# Pseudobulk correlation, i.e. average TFA over cells of the same cell type
###
pseudobulk_mean <- function(mat, celltypes, cellnames) {
    sapply(celltypes, function(ct) {
        cells_ct <- intersect(
            cellnames[[ct]],
            colnames(mat)
        )
        if (length(cells_ct) == 0) {
            return(rep(NA_real_, nrow(mat)))
        }
        rowMeans(
            mat[, cells_ct, drop = FALSE],
            na.rm = TRUE
        )
    })
}

chrom_pseudobulk <- pseudobulk_mean(
    as.matrix(cell_reordered_filtered_chrom_assay),
    celltypes,
    cellnames
)


rownames(chrom_pseudobulk) <- all_reg
system('mkdir -p TFA_output/analysis/TFA_corr/pseudobulk')

results_corr_psb <- mclapply(algos, function(a){
    cur_P <- as.matrix(get_concate_P(a))
    rownames(cur_P) <- all_reg
    colnames(cur_P) <- Ps_cellnames
    TFA_pseudobulk <- pseudobulk_mean(
        cur_P,
        celltypes,
        cellnames
    )
    pearson_correlations <- sapply(
        seq_len(nrow(TFA_pseudobulk)),
        function(i) {
            cor(
                TFA_pseudobulk[i, ],
                chrom_pseudobulk[i, ],
                method = "pearson"
            )
        }
    )
    spearman_correlations <- sapply(
        seq_len(nrow(TFA_pseudobulk)),
        function(i) {
            cor(
                TFA_pseudobulk[i, ],
                chrom_pseudobulk[i, ],
                method = "spearman"
            )
        }
    )
    write.table(
        pearson_correlations,
        paste0(
            'TFA_output/analysis/TFA_corr/pseudobulk/Pearson_',
            a, '.txt'
        ),
        quote=FALSE,
        col.names=FALSE,
        row.names=FALSE
    )
    write.table(
        spearman_correlations,
        paste0(
            'TFA_output/analysis/TFA_corr/pseudobulk/Spearman_',
            a, '.txt'
        ),
        quote=FALSE,
        col.names=FALSE,
        row.names=FALSE
    )
    return(list(
        pearson = pearson_correlations,
        spearman = spearman_correlations
    ))
}, mc.cores=max_cores)

pearson_lst <- setNames(lapply(seq_along(algos), function(i) results_corr_psb[[i]]$pearson), algos)
spearman_lst <- setNames(lapply(seq_along(algos), function(i) results_corr_psb[[i]]$spearman), algos)

pearson_tab <- data.frame(
  algo=unlist(lapply(algos, function(a) rep(a, length(all_reg)))),
  pearson=unlist(pearson_lst)
)
pearson_tab <- formal_algo_names(pearson_tab, "algo")

spearman_tab <- data.frame(
  algo=unlist(lapply(algos, function(a) rep(a, length(all_reg)))),
  spearman=unlist(spearman_lst)
)
spearman_tab <- formal_algo_names(spearman_tab, "algo")

p_pearson_box <- ggplot(pearson_tab, aes(x=algo, y=pearson)) +
  geom_boxplot() +
  ggtitle("Pearson correlation coefficients across pseudobulk cell types") +
  theme_minimal() +
  theme( axis.text.x=element_text(angle=45, hjust=1, size=14))

ggsave("TFA_output/analysis/TFA_corr/pseudobulk/Pearson.pdf", p_pearson_box)

p_spearman_box <- ggplot(spearman_tab, aes(x=algo, y=spearman)) +
  geom_boxplot() +
  ggtitle("Spearman correlation coefficients across pseudobulk cell types") +
  theme_minimal() +
  theme( axis.text.x=element_text(angle=45, hjust=1, size=14))

ggsave("TFA_output/analysis/TFA_corr/pseudobulk/Spearman.pdf", p_spearman_box)

pearson_group <- aggregate(pearson ~ algo, pearson_tab, function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE)))
pearson_group <- data.frame(algo=pearson_group$algo, mean=pearson_group$pearson[, "mean"], sd=pearson_group$pearson[, "sd"])

spearman_group <- aggregate(spearman ~ algo, spearman_tab, function(x) c(mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE)))
spearman_group <- data.frame(algo=spearman_group$algo, mean=spearman_group$spearman[, "mean"], sd=spearman_group$spearman[, "sd"])

excluded_algo <- c("HALS-p", "IBPG-p")
pearson_group <- pearson_group[!pearson_group$algo %in% excluded_algo, ]
spearman_group <- spearman_group[!spearman_group$algo %in% excluded_algo, ]

p1 <- ggplot(pearson_group, aes(algo, mean)) +
  geom_col(fill="skyblue", color="black", width=0.7) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, color="red") +
  geom_text(aes(label=round(mean, 4)), vjust=-0.5, size=5) +
  labs(title="Global Pearson CC across pseudobulk cell types", x="", y="") +
  theme_minimal() +
  theme( axis.text.x=element_text(angle=45, hjust=1, size=14), axis.text.y=element_text(size=12), plot.title=element_text(size=16))

ggsave("TFA_output/analysis/TFA_corr/pseudobulk/mean_Pearson.pdf", p1)

p2 <- ggplot(spearman_group, aes(algo, mean)) +
  geom_col(fill="skyblue", color="black", width=0.7) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, color="red") +
  geom_text(aes(label=round(mean, 4)), vjust=-0.5, size=5) +
  labs(title="Global Spearman CC across pseudobulk cell types", x="", y="") +
  theme_minimal() +
  theme( axis.text.x=element_text(angle=45, hjust=1, size=14), axis.text.y=element_text(size=12), plot.title=element_text(size=16))

ggsave("TFA_output/analysis/TFA_corr/pseudobulk/mean_Spearman.pdf", p2)

###
# 3. plotting TFA ("TFAct": default 500, 100)
###
library(Seurat)

pbmc <- readRDS("TIGER_zenodo/data/PBMC/pbmc_final.rds")

ct_algo <- c("ibpg_2fct_500_0", "ibpg_2fct", "ibpg_2fct_500_1000")
ct_concatenated_P <- mclapply(ct_algo, function(a){
    cur_P <- get_concate_P(a)
    colnames(cur_P) <- colnames(cell_reordered_filtered_chrom_assay)
    rownames(cur_P) <- rownames(cell_reordered_filtered_chrom_assay)
    # rearrange to make the columns (cells) follow the same order of pbmc 
    cur_P <- cur_P[,colnames(pbmc)]
    return(cur_P)
    }, mc.cores=20)

ct_algo_name_corrected <- c("ibpg_2fct_500_0", "ibpg_2fct_500_100", "ibpg_2fct_500_1000")
names(ct_concatenated_P) <- ct_algo_name_corrected


for(a in ct_algo_name_corrected){
    pbmc2 <- pbmc
    alg_name <- a
    cur_assay = CreateAssayObject(data=as.matrix(ct_concatenated_P[[alg_name]]))
    pbmc2[[alg_name]] = cur_assay

    cur_out <- paste0('plots/TFA/', alg_name)
    system(paste0('mkdir -p plots/TFA/', alg_name))

    mclapply(rownames(pbmc2[[alg_name]]), function(f){
        DefaultAssay(pbmc2) = alg_name
        FeaturePlot(pbmc2, features = f, reduction='wnn.umap')
        ggsave(paste0(cur_out, '/', f, '.png'), width = 10, height = 8)
    }, mc.cores=50)
}

# see if this is different
DefaultAssay(pbmc2) <- 'RNA'
system('mkdir -p plots/TFA/RNA')
mclapply(all_reg, function(f){
        FeaturePlot(pbmc2, features = f, reduction='wnn.umap')
        ggsave(paste0('plots/TFA/RNA/', f, '.png'), width = 10, height = 8)
    }, mc.cores=50)

# DE and DA TFs
system(paste0('mkdir -p plots/DE/'))
DE_gene_lst <- list()
for(alg_name in ct_algo_name_corrected){
    cur_assay = CreateAssayObject(data=as.matrix(ct_concatenated_P[[alg_name]]))
    pbmc2[[alg_name]] = cur_assay
    DefaultAssay(pbmc2) = alg_name
    all_markers <- FindAllMarkers(pbmc2, group.by="celltype")
    DE_gene_lst[[length(DE_gene_lst)+1]] = all_markers
    write.table(all_markers, paste0('plots/DE/', alg_name, '.txt'), quote=FALSE, col.names = TRUE, row.names=FALSE,sep='\t')
}

names(DE_gene_lst) <- ct_algo_name_corrected

DefaultAssay(pbmc2) = "RNA"
all_markers <- FindAllMarkers(pbmc2, features = all_reg, group.by="celltype")
write.table(all_markers, paste0('plots/DE/RNA.txt'), quote=FALSE, col.names = TRUE, row.names=FALSE,sep="\t")

# # Use the following to find a list of TF significantly upregulated in cell type HSPC
# all_markers[all_markers$cluster=='HSPC' & all_markers$p_val_adj<0.1 & all_markers$avg_log2FC > 0,][['gene']]

# # Use the following to find a list of TF significantly more active in cell type HSPC
# cur_markers = DE_gene_lst[['ibpg_2fct_500_100']]
# cur_markers[cur_markers$cluster=='HSPC'&cur_markers$p_val_adj<0.1 & cur_markers$avg_log2FC > 0,][['gene']]

cur_markers = DE_gene_lst[['ibpg_2fct_500_100']]
### Make a list of upregulated / more active TFs, or downregulated / less active TFs
lapply(names(DE_gene_lst), function(n){system(paste0('mkdir -p plots/DE/up_regulated/', n))})
system(paste0('mkdir -p plots/DE/up_regulated/RNA'))

lapply(names(DE_gene_lst), function(n){system(paste0('mkdir -p plots/DE/down_regulated/', n))})
system(paste0('mkdir -p plots/DE/down_regulated/RNA'))

cell_clust = unique(cur_markers$cluster)
for(c in cell_clust){
    c_sub = gsub(' ','_',c)
    for(n in names(DE_gene_lst)){
        cur_markers = DE_gene_lst[[n]]
        gene <- cur_markers[cur_markers$cluster==c&cur_markers$p_val_adj<0.1 & cur_markers$avg_log2FC > 0,][['gene']]
        write.table(gene, paste0('plots/DE/up_regulated/', n,'/',c_sub,'.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
        gene <- cur_markers[cur_markers$cluster==c&cur_markers$p_val_adj<0.1 & cur_markers$avg_log2FC < 0,][['gene']]
        write.table(gene, paste0('plots/DE/down_regulated/', n,'/',c_sub,'.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
    }
    gene <- all_markers[all_markers$cluster==c&all_markers$p_val_adj<0.1 & all_markers$avg_log2FC > 0,][['gene']]
    write.table(gene, paste0('plots/DE/up_regulated/RNA/',c_sub,'.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
    gene <- all_markers[all_markers$cluster==c&all_markers$p_val_adj<0.1 & all_markers$avg_log2FC < 0,][['gene']]
    write.table(gene, paste0('plots/DE/down_regulated/RNA/',c_sub,'.txt'), quote=FALSE, col.names = FALSE, row.names=FALSE)
}

DATF_outdir="TFA_output/analysis/DATF"
system('mkdir -p TFA_output/analysis/DATF')

for (cur_name in names(DE_gene_lst)){
    DA_TF_and_or_not_DE <- mclapply(cell_clust, function(c){
        cur_markers = DE_gene_lst[[cur_name]]
        TF_A <- cur_markers[cur_markers$cluster==c&cur_markers$p_val_adj<0.1 & cur_markers$avg_log2FC > 0,][['gene']]
        TF_up <- all_markers[all_markers$cluster==c&all_markers$p_val_adj<0.1 & all_markers$avg_log2FC > 0,][['gene']]
        return(list(not_up = setdiff(TF_A, TF_up), up = intersect(TF_A, TF_up)))
    }, mc.cores=max_cores)

    DA_TF_not_DE <- list()
    DA_TF_and_DE <- list()

    for (i in seq_along(cell_clust)) {
        DA_TF_not_DE[[i]] <- DA_TF_and_or_not_DE[[i]]$not_up
        DA_TF_and_DE[[i]] <- DA_TF_and_or_not_DE[[i]]$up
    }

    names(DA_TF_not_DE) <- cell_clust
    names(DA_TF_and_DE) <- cell_clust

    # for DA_TF_not_DE TFs, we further filter for TFs that at least have some expression in the corresponding cell type

    clust_expr <- AverageExpression(pbmc, assays = 'RNA', layer = "counts", group.by='celltype')
    clust_expr <- (clust_expr[[1]])[all_reg,]

    # average expression >= 0.1 (assuming high dropout)
    colnames(clust_expr) <- gsub('-', '_', colnames(clust_expr))

    DA_TF_not_DE_filtered <- lapply(names(DA_TF_not_DE), function(c){
        cur_expr <- sapply(DA_TF_not_DE[[c]], function(tf){clust_expr[tf,c]})
        return(names(cur_expr[cur_expr>=0.1]))
        })

    names(DA_TF_not_DE_filtered) <- names(DA_TF_not_DE)

    DATF_df_up <- data.frame(cell=names(DA_TF_and_DE), DATF=sapply(DA_TF_and_DE, function(x){paste(x,collapse=', ')}))
    write.table(DATF_df_up, paste0(DATF_outdir,'/',cur_name,"_DATF_df_up.txt"), quote=FALSE, sep = "\t",row.names = FALSE)

    DATF_df <- data.frame(cell=names(DA_TF_not_DE_filtered), DATF=sapply(DA_TF_not_DE_filtered, function(x){paste(x,collapse=', ')}))
    write.table(DATF_df, paste0(DATF_outdir,'/',cur_name,"_DATF_df.txt"), quote=FALSE, sep = "\t",row.names = FALSE)

}

## concensus DATF
DATF_df_up_500_100 <- read.table(paste0(DATF_outdir,"/ibpg_2fct_500_100_DATF_df_up.txt"), sep = "\t", header=TRUE)
DATF_df_500_100 <- read.table(paste0(DATF_outdir,"/ibpg_2fct_500_100_DATF_df.txt"), sep = "\t",header=TRUE)
DATF_df_up_500_1000 <- read.table(paste0(DATF_outdir,"/ibpg_2fct_500_1000_DATF_df_up.txt"), sep = "\t", header=TRUE)
DATF_df_500_1000 <- read.table(paste0(DATF_outdir,"/ibpg_2fct_500_1000_DATF_df.txt"), sep = "\t",header=TRUE)

DATF_df_intersected <- lapply(seq(1,length(cell_clust)), function(x){intersect(strsplit(DATF_df_500_100$DATF, ', ')[[x]], strsplit(DATF_df_500_1000$DATF, ', ')[[x]])})
DATF_df_up_intersected <- lapply(seq(1,length(cell_clust)), function(x){intersect(strsplit(DATF_df_up_500_100$DATF, ', ')[[x]], strsplit(DATF_df_up_500_1000$DATF, ', ')[[x]])})
names(DATF_df_intersected) <- cell_clust
names(DATF_df_up_intersected) <- cell_clust
DATF_df_intersected_df <- data.frame(cell=names(DATF_df_intersected), DATF=sapply(DATF_df_intersected, function(x){paste(x,collapse=', ')}))
DATF_df_up_intersected_df <- data.frame(cell=names(DATF_df_up_intersected), DATF=sapply(DATF_df_up_intersected, function(x){paste(x,collapse=', ')}))
write.table(DATF_df_intersected_df, paste0(DATF_outdir,"/100_1000_intersected_DATF_df.txt"), quote=FALSE, sep = "\t",row.names = FALSE)
write.table(DATF_df_up_intersected_df, paste0(DATF_outdir,"/100_1000_intersected_DATF_df_up.txt"), quote=FALSE, sep = "\t",row.names = FALSE)


DE_TFs <- lapply(cell_clust, function(c){all_markers[all_markers$cluster==c&all_markers$p_val_adj<0.1 & all_markers$avg_log2FC > 0,][['gene']]})
names(DE_TFs) <- cell_clust
