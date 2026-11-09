#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data

# use rep 0
library(ggplot2)
library(parallel)

outdir="TFA_out"

datasets=c("processed_FBS_reprogramming", "processed_Buenrostro_Hematopoiesis", "processed_Ranzoni_fetal_hematopoiesis_fine_tree")
selected_cluster=c("cluster4", "cluster7", "HSC-MPPs") # largest clusters
names(selected_cluster) <- datasets
algo_E_100=c("ibpg_2fct", "ibpg_2fct_2000_100") # former: 500, 100
algo_E_1000=c("ibpg_2fct_100_1000", "ibpg_2fct_2000_1000")
algos=c(algo_E_100,algo_E_1000)
rep_cnt=20

formal_algo_names <- function(dat, column_name){
    dat[[column_name]][dat[[column_name]] == "ibpg"] <- "IBPG-p"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct"] <- "TFAct-F_500_100"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_1000"] <- "TFAct-F_500_1000"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_500_0"] <- "TFAct-F_500_0"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_2000_100"] <- "TFAct-F_2000_100"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_2000_1000"] <- "TFAct-F_2000_1000"
    dat[[column_name]][dat[[column_name]] == "ibpg_2fct_100_1000"] <- "TFAct-F_100_1000"
    return(dat)
}

count_nonzero <- function(f_name){
    tab <- read.table(f_name, header=FALSE)
    return(sum(tab!=0))
}
cur_A_lst <- list()
for(a in algos){
    for(i in 1:length(datasets)){
        cur_counts <- unlist(mclapply(0:(rep_cnt-1), function(r){count_nonzero(paste0("TFA_output/", a, "/", datasets[[i]], "/", selected_cluster[[i]],"/rep_", r, "/A.txt"))}, mc.cores=30))
        cur_A_lst[[length(cur_A_lst)+1]] <- data.frame(algo=a, ds=datasets[[i]], num_nonzero=cur_counts)
    }
}

A_tab <- do.call(rbind, cur_A_lst)

A_tab <- formal_algo_names(A_tab, "algo")
E_100=c("TFAct-F_500_100", "TFAct-F_2000_100") # former: 500, 100
E_1000=c("TFAct-F_100_1000", "TFAct-F_2000_1000")

system("mkdir -p TFA_output/analysis/sparsity")

for(d in datasets){
    cur = A_tab[A_tab$ds==d,]
    ggplot(cur, aes(x = factor(algo), y = num_nonzero)) +
    geom_boxplot() + 
    labs(
        title = "number of nonzero entries",
        x = "",
        y = "",
        color = "algo"
    ) +
    theme_minimal() + theme(text=element_text(size=15))
    ggsave(paste0("TFA_output/analysis/sparsity/",d,'_',selected_cluster[[d]],'.png'), width = 8, height = 6, dpi = 300)
}

# TFA_output/ibpg_2fct/processed_Buenrostro_Hematopoiesis/cluster7/rep_0/A.txt
