#pwd=/mnt/dv/wid/projects6/Roy-singlecell3/yuda_work/TFA_estimation/real_dat_validation/scMTNI_data
library(ggplot2)

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

PTF <- read.table("predictable_tf_nums.txt", header=FALSE, check.names=FALSE, sep='\t')

colnames(PTF) <- c('algo', 'ds', 'c', 'gold', 'num_predictable_TFs')

PTF <- formal_algo_names(PTF, "algo")

system('mkdir -p TFA_output/analysis/predictable_TFs')
# FBS
FBS <- PTF[PTF$ds=="FBS",]
for(gs in unique(FBS$gold)){
    cur_df <- FBS[FBS$gold==gs,]
    ggplot(cur_df, aes(y=num_predictable_TFs, x=algo)) + geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), text = element_text(size = 15))
    ggsave(paste0('TFA_output/analysis/predictable_TFs/FBS_',gs,'.png'))
}

# BH
BH <- PTF[PTF$ds=="BH",]
for(gs in unique(BH$gold)){
    cur_df <- BH[BH$gold==gs,]
    ggplot(cur_df, aes(y=num_predictable_TFs, x=algo)) + geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), text = element_text(size = 15))
    ggsave(paste0('TFA_output/analysis/predictable_TFs/BH_',gs,'.png'))
}

# RFH
RFH <- PTF[PTF$ds=="RFH",]
for(gs in unique(RFH$gold)){
    cur_df <- RFH[RFH$gold==gs,]
    ggplot(cur_df, aes(y=num_predictable_TFs, x=algo)) + geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), text = element_text(size = 15))
    ggsave(paste0('TFA_output/analysis/predictable_TFs/RFH_',gs,'.png'))
}