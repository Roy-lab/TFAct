library(ggplot2)
library("scmamp")

datasets <- c('processed_FBS_reprogramming', 'processed_Buenrostro_Hematopoiesis', 'processed_Ranzoni_fetal_hematopoiesis_fine_tree')
ignored_algo <- c('IBPG-p', 'HALS-p')

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

for (ds in datasets){

    Fscore_file <- paste0("TFA_output/",ds,"/analysis/Fscores.txt")

    F_tab <-read.table(Fscore_file, header=TRUE, sep = "\t", check.names=FALSE)

    F_tab$algo <- gsub("^.*/([^/]*)$","\\1", F_tab$algo)

    F_tab <- formal_algo_names(F_tab, "algo")

    F_tab$Fscore_top_edge <- as.factor(F_tab$Fscore_top_edge)

	F_tab <- F_tab[!(F_tab$algo %in% ignored_algo), ]

	if (ds != "processed_Ranzoni_fetal_hematopoiesis_fine_tree"){
		F_tab$cell <- gsub("^cluster", "C", F_tab$cell)
	}

    # for(E in unique(F_tab$gold_standard)){
	# 	cur_F_tab <- F_tab[F_tab$gold_standard == E,]
	# 	ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_gt, x=Fscore_top_edge)) +
	# 		geom_boxplot(position=position_dodge(width = 0.8)) +
	# 	theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))
	# 	ggsave(paste0("TFA_output/",ds,"/analysis/Fscores_ground_true_", gsub("^E", "", E), ".png"), width = 15, height = 10)

	# 	ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_prior, x=Fscore_top_edge)) +
	# 		geom_boxplot(position=position_dodge(width = 0.8)) +
	# 	theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))+
	# 	scale_x_discrete(expand = c(0, -0.8))
	# 	ggsave(paste0("TFA_output/",ds,"/analysis/Fscores_prior_", gsub("^E", "", E), ".png"), width = 15, height = 10)
    # }

	# top 1000
	system(paste0('mkdir -p ', "TFA_output/",ds,"/analysis/gt"))
	system(paste0('mkdir -p ', "TFA_output/",ds,"/analysis/prior"))
	for(E in unique(F_tab$Fscore_top_edge)){	
		for(gs in unique(F_tab$gold_standard)){
			cur_F_tab <- F_tab[F_tab$Fscore_top_edge == E & F_tab$gold_standard==gs,]
			cur_F_tab$fscore_w_gt <- as.numeric(cur_F_tab$fscore_w_gt)

			midpoint <- (min(cur_F_tab$fscore_w_gt, na.rm = TRUE) + max(cur_F_tab$fscore_w_gt, na.rm = TRUE)) / 2
			cur_F_tab$text_color <- ifelse(cur_F_tab$fscore_w_gt > midpoint, "white", "black")

			ggplot(cur_F_tab, aes(x = cell, y = algo, fill = fscore_w_gt)) +
				geom_tile(color = "white") +
				scale_fill_gradient(low = "white",high = "blue",limits = range(cur_F_tab$fscore_w_gt, na.rm = TRUE)) +
				geom_text(aes(label = round(fscore_w_gt, 3), color = text_color), size = 6) +
				scale_color_identity() +
				labs(title = "", x = "", y = "", fill = "F-score") +
				theme_minimal() +
				theme(
					axis.text.x = element_text(angle = 90, hjust = 1),
					text = element_text(size = 20)
				)

			# ggplot(cur_F_tab, aes(y=fscore_w_gt, x=algo)) +
			# 	geom_boxplot(position=position_dodge(width = 0.3), width = 0.5) +
			# theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size = 20)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
			if (ds == "processed_Ranzoni_fetal_hematopoiesis_fine_tree") {
				ggsave(paste0("TFA_output/",ds,"/analysis/gt/", gs,"_top", E, ".pdf"), width = 12, height = 7)
			} else {
				ggsave(paste0("TFA_output/",ds,"/analysis/gt/", gs,"_top", E, ".pdf"), width = 10, height = 6.5)
			}

			cur_F_tab_df <- reshape(cur_F_tab[,c("cell", "algo", "fscore_w_gt")], idvar="cell", timevar="algo", direction="wide")
			colnames(cur_F_tab_df) <- gsub("fscore_w_gt.", "", colnames(cur_F_tab_df))
			cur_F_tab_df[,"cell"] <- NULL

			pdf(paste0("TFA_output/",ds,"/analysis/gt/critical_difference_", gs,"_top", E, ".pdf"), width = 10, height = 6)
			plotCD(cur_F_tab_df, alpha=0.05, cex=0.5, lwd = 2)
			dev.off()

			cur_F_tab$fscore_w_prior <- as.numeric(cur_F_tab$fscore_w_prior)
			midpoint <- (min(cur_F_tab$fscore_w_prior, na.rm = TRUE) + max(cur_F_tab$fscore_w_prior, na.rm = TRUE)) / 2

			cur_F_tab$text_color <- ifelse(cur_F_tab$fscore_w_prior > midpoint, "white", "black")

			ggplot(cur_F_tab, aes(x = cell, y = algo, fill = fscore_w_prior)) +
				geom_tile(color = "white") +
				scale_fill_gradient(low = "white",high = "blue",limits = range(cur_F_tab$fscore_w_prior, na.rm = TRUE)) +
				geom_text(aes(label = round(fscore_w_prior, 3), color = text_color), size = 6) +
				scale_color_identity() +
				labs(title = "", x = "", y = "", fill = "F-score") +
				theme_minimal() +
				theme(
					axis.text.x = element_text(angle = 90, hjust = 1),
					text = element_text(size = 20)
				)

			# ggplot(cur_F_tab, aes(y=fscore_w_prior, x=algo)) +
			# 	geom_boxplot(position=position_dodge(width = 0.8), width = 0.5) +
			# theme(axis.text.x = element_text(angle = 90, hjust = 1),text = element_text(size = 20))+
			# scale_x_discrete(expand = c(0, -0.8)) + theme(axis.title.x=element_blank())
			if (ds == "processed_Ranzoni_fetal_hematopoiesis_fine_tree") {
				ggsave(paste0("TFA_output/",ds,"/analysis/prior/", gs,"_top",  E, ".pdf"), width = 12, height = 7)
			} else {
				ggsave(paste0("TFA_output/",ds,"/analysis/prior/", gs,"_top",  E, ".pdf"), width = 10, height = 6.5)
			}
		}
		# cur_F_tab <- F_tab[F_tab$Fscore_top_edge == E,]
		# ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_gt, x=gold_standard)) +
		# 	geom_boxplot(position=position_dodge(width = 0.8)) +
		# theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))
		# ggsave(paste0("TFA_output/",ds,"/analysis/Fscores_ground_true_top", E, ".png"), width = 15, height = 10)

		# ggplot(cur_F_tab, aes(fill=algo, y=fscore_w_prior, x=gold_standard)) +
		# 	geom_boxplot(position=position_dodge(width = 0.8)) +
		# theme(axis.text.x = element_text(angle = 45, hjust = 1),text = element_text(size = 20))+
		# scale_x_discrete(expand = c(0, -0.8))
		# ggsave(paste0("TFA_output/",ds,"/analysis/Fscores_prior_",  E, ".png"), width = 15, height = 10)
	}

}


