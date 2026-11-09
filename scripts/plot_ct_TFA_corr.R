library(ggplot2)
library(parallel)


formal_algo_names <- function(dat, column_name){
  for(i in c(0, 10, 100, 1000)){
    for(j in c(0, 10, 100, 1000)){
      dat[[column_name]][dat[[column_name]] == paste("ibpg_2fct", i, j,sep="_")] <- paste("F", i, j,sep="_")
      dat[[column_name]][dat[[column_name]] == paste("kl_nmf_2fct", i, j,sep="_")] <- paste("KL", i, j,sep="_")
    }
  }   
  return(dat)
}

algos=c('ibpg_2fct', 'kl_nmf_2fct')

#algos=c('baselines/hals_prior', 'FastNCA', 'NINCA_quad', 'pmf_grn', 'regularized_nca', 'BITFAM', 'TIGER', 'NetREX_CF', 'supirfactor_shallow', 'supirfactor_hierarchical', 'ibpg_2f')
Edge_adding_probs=c(0.2, 0.4) # edge retaining probability
Ldaks=c(10, 20)
Prefixes=c('Spom', 'Anc2', 'Anc5', 'Calb', 'Cgla', 'right', 'Sbay')
Exp_mat_w_dropouts=c('E', 'E_0.6')
# Exp_mat_w_dropouts=c('E', 'E_0.3', 'E_0.6')
TFA_sim_reps=c('rep0')
l1_choices = c(0, 10, 100, 1000)
lineage_reg_choices = c(0, 10, 100, 1000)
l_combo = expand.grid(l1_choices, lineage_reg_choices)
l_combo <- paste(l_combo$Var1, l_combo$Var2, sep='_')

result<- expand.grid(Edge_adding_probs, Ldaks, Prefixes, Exp_mat_w_dropouts, TFA_sim_reps)
list_of_vectors <- split(result, seq(nrow(result)))

plot_TFA <- function(edge_adding_prob, ldak, prefix, exp_mat, sim_rep){

    combined_df <-data.frame()
    for(i in 1:length(algos)){
        for (j in 1:length(l_combo)){
            cur_f <- paste0("TFA_output/ct_algos/",algos[i],"/",l_combo[j], "/", edge_adding_prob,"/k",ldak,"/",prefix,"/rep0/",exp_mat,"_stability/stats_1000/tf_specifics.txt") # note that stats_1000, stats_2000's tf_specifics are the same as the number just affects fscores.
            cur_df <-read.table(cur_f, header=TRUE, sep = "\t", check.names=FALSE)

            cur_df$algo <- paste(gsub("^.*/([^/]*)$","\\1",algos[i]), l_combo[j], sep='_')

            combined_df <- rbind(combined_df, cur_df)
        }
    }

    combined_df <- formal_algo_names(combined_df, "algo")

    suffix=paste0(edge_adding_prob, "_k", ldak,"_",prefix,"_",sim_rep,"_",exp_mat,"_1000.pdf")
    # suffix=paste0(edge_adding_prob, "_k", ldak,"_",prefix,"_",sim_rep,"_",exp_mat,"_1000.png")

    ggplot(combined_df, aes(y=Pearson, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 30), text = element_text(size = 15))

    ggsave(paste0("TFA_output/analysis/ct/tfa_corr/Pearson_",suffix))
    ggplot(combined_df, aes(y=Spearman, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + theme(axis.title.x=element_blank(), axis.title.y=element_blank()) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 30), text = element_text(size = 15))

    ggsave(paste0("TFA_output/analysis/ct/tfa_corr/Spearman_",suffix))

    return(combined_df[,c('algo', 'Pearson', 'Spearman')])
}

plot_all_TFA_by_group <- function(tab, suffix){
    tab <- formal_algo_names(tab, "algo")
    ggplot(tab, aes(y=Pearson, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  axis.text.y = element_text(size = 30), text = element_text(size = 13)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

    ggsave(paste0("TFA_output/analysis/ct/tfa_corr/Pearson_",suffix))
    ggplot(tab, aes(y=Spearman, x=algo)) + 
	geom_boxplot(position=position_dodge(width = 0.8)) + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  axis.text.y = element_text(size = 30), text = element_text(size = 13))+ 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())

    ggsave(paste0("TFA_output/analysis/ct/tfa_corr/Spearman_",suffix))
}

outdir='TFA_output/analysis/ct/tfa_corr'
system(paste0('mkdir -p ', outdir))

stat_lst<-mclapply(list_of_vectors, function(x){plot_TFA(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]])}, mc.cores = 20)

# covariates that may affect TFA
cov_classes <- lapply(list_of_vectors, function(x){paste(x$Var1, x$Var2, x$Var4, sep=",")}) # Var1 is 0.2/0.4, Var2 is 20/10, Var4 is expression, potentially with dropouts

unique_classes <- unique(cov_classes)

mclapply(c(1:length(unique_classes)), function(x){plot_all_TFA_by_group(do.call(rbind, stat_lst[cov_classes %in% unique_classes[[x]]]), paste0(gsub(",", "_", unique_classes[[x]]), ".pdf"))})
# mclapply(c(1:length(unique_classes)), function(x){plot_all_TFA_by_group(do.call(rbind, stat_lst[cov_classes %in% unique_classes[[x]]]), paste0(gsub(",", "_", unique_classes[[x]]), ".png"))})

###
# Expanded hyperparameter search. By ChatGPT
###
################################################################################
# Expanded x/y analysis
################################################################################

# Change this one line if you want a different font in Illustrator
FONT_FAMILY <- "Arial"

expanded_outdir <- "TFA_output/analysis/ct/tfa_corr_expanded"
dir.create(expanded_outdir, recursive=TRUE, showWarnings=FALSE)

xy_vals <- c(
    "0",
    "10",
    "100",
    "1000",
    "10000",
    "100000",
    "1000000",
    "10000000"
)

# Shorter labels make the figures substantially less crowded
xy_labels <- c(
    "0"        = "0",
    "10"       = "10",
    "100"      = "100",
    "1000"     = "1k",
    "10000"    = "10k",
    "100000"   = "100k",
    "1000000"  = "1M",
    "10000000" = "10M"
)


################################################################################
# Read expanded results
################################################################################

read_expanded_TFA <- function(edge_adding_prob, ldak, prefix, exp_mat, sim_rep){

    dfs <- list()
    counter <- 1

    for(algo in algos){
        for(x in xy_vals){
            for(y in xy_vals){

                cur_f <- paste0(
                    "TFA_output/ct_algos/",
                    algo, "/",
                    x, "_", y, "/",
                    edge_adding_prob,
                    "/k", ldak, "/",
                    prefix, "/",
                    sim_rep, "/",
                    exp_mat,
                    "_stability/stats_1000/tf_specifics.txt"
                )

                # skip missing runs rather than stopping the whole script
                if(!file.exists(cur_f)){
                    next
                }

                cur_df <- read.table(
                    cur_f,
                    header=TRUE,
                    sep="\t",
                    check.names=FALSE
                )

                cur_df$algo <- ifelse(
                    algo == "ibpg_2fct",
                    "F",
                    "KL"
                )

                cur_df$x <- x
                cur_df$y <- y
                cur_df$prefix <- prefix
                cur_df$edge_adding_prob <- edge_adding_prob
                cur_df$ldak <- ldak
                cur_df$exp_mat <- exp_mat

                dfs[[counter]] <- cur_df
                counter <- counter + 1
            }
        }
    }

    if(length(dfs) == 0){
        return(NULL)
    }

    do.call(rbind, dfs)
}

Edge_adding_probs_expanded <- c(0.4)
Ldaks_expanded <- c(20)

result_expanded <- expand.grid(
    Edge_adding_probs_expanded,
    Ldaks_expanded,
    Prefixes,
    Exp_mat_w_dropouts,
    TFA_sim_reps
)

list_of_vectors_expanded <- split(
    result_expanded,
    seq(nrow(result_expanded))
)

expanded_stat_lst <- mclapply(
    list_of_vectors_expanded,
    function(z){
        read_expanded_TFA(
            z[[1]],
            z[[2]],
            z[[3]],
            z[[4]],
            z[[5]]
        )
    },
    mc.cores=20
)

expanded_stat_lst <- expanded_stat_lst[
    !sapply(expanded_stat_lst, is.null)
]

expanded_tab <- do.call(rbind, expanded_stat_lst)


################################################################################
# Publication-quality theme
################################################################################

expanded_theme <- theme_classic(
    base_family=FONT_FAMILY,
    base_size=18
) +
    theme(
        axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=20),
        strip.text=element_text(size=18),
        strip.background=element_blank(),
        panel.spacing=grid::unit(0.6, "lines"),
        plot.margin=margin(8, 8, 8, 8)
    )


################################################################################
# Boxplots
#
# One panel for each algorithm x x-value.
# Within each panel, the x-axis is the changing y value.
################################################################################

plot_expanded_boxplot <- function(dat, metric, outfile){

    dat$x <- factor(dat$x, levels=xy_vals)
    dat$y <- factor(dat$y, levels=xy_vals)

    p <- ggplot(
        dat,
        aes(
            x=y,
            y=.data[[metric]]
        )
    ) +
        geom_boxplot(
            width=0.70,
            outlier.size=0.7,
            linewidth=0.35,
            na.rm=TRUE
        ) +
        facet_grid(
            algo ~ x,
            drop=FALSE,
            labeller=labeller(
                x=function(z) paste0("lambda[C] == '", xy_labels[z], "'"),
                .default=label_parsed
            )
        ) +
        scale_x_discrete(
            labels=xy_labels,
            drop=FALSE
        ) +
        labs(
            x=expression(lambda[E]),
            y=metric
        ) +
        expanded_theme +
        theme(
            axis.text.x=element_text(
                angle=45,
                hjust=1,
                vjust=1,
                size=14
            ),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=19),
            strip.text=element_text(size=17),
            panel.spacing.x=grid::unit(0.3, "lines")
        )

    ggsave(
        outfile,
        p,
        width=14,
        height=7,
        device=cairo_pdf
    )
}


################################################################################
# Heatmaps
#
# Median is calculated from exactly the values entering the boxplots.
# Two subplots: F and KL.
################################################################################

plot_expanded_heatmap <- function(dat, metric, outfile){

    med_tab <- aggregate(
        dat[[metric]],
        by=list(
            algo=dat$algo,
            x=dat$x,
            y=dat$y
        ),
        FUN=function(z){
            z <- z[is.finite(z)]
            if(length(z) == 0) NA_real_ else median(z)
        }
    )

    colnames(med_tab)[4] <- "value"

    # Make sure every x-y combination exists, even if all values were NA
    full_grid <- expand.grid(
        algo=unique(dat$algo),
        x=xy_vals,
        y=xy_vals
    )

    med_tab <- merge(
        full_grid,
        med_tab,
        by=c("algo", "x", "y"),
        all.x=TRUE
    )

    med_tab$x <- factor(
        med_tab$x,
        levels=xy_vals,
        labels=xy_labels
    )

    med_tab$y <- factor(
        med_tab$y,
        levels=rev(xy_vals),
        labels=rev(xy_labels)
    )

    p <- ggplot(
        med_tab,
        aes(
            x=x,
            y=y,
            fill=value
        )
    ) +
        geom_tile(
            colour="white",
            linewidth=0.35
        ) +
        geom_text(
            aes(label=ifelse(is.na(value), "NA", sprintf("%.2f", value))),
            family=FONT_FAMILY,
            size=5.5
        ) +
        facet_wrap(
            ~algo,
            nrow=1
        ) +
        scale_fill_gradient2(
            low="blue",
            mid="white",
            high="red",
            midpoint=0,
            limits=c(-1, 1),
            na.value="grey90",
            name=paste0("Median\n", metric)
        ) +
        labs(
            x=expression(lambda[C]),
            y=expression(lambda[E])
        ) +
        coord_equal() +
        expanded_theme +
        theme(
            axis.text.x=element_text(
                angle=45,
                hjust=1
            ),
            legend.title=element_text(size=11),
            legend.text=element_text(size=10)
        )

    ggsave(
        outfile,
        p,
        width=12,
        height=6,
        device=cairo_pdf
    )
    
}

plot_expanded_NA_heatmap <- function(dat, metric, outfile){

    na_tab <- aggregate(
        dat[[metric]],
        by=list(
            algo=dat$algo,
            x=dat$x,
            y=dat$y
        ),
        FUN=function(z){
            sum(!is.finite(z))
        }
    )

    colnames(na_tab)[4] <- "n_na"

    # Ensure all x/y combinations are represented
    full_grid <- expand.grid(
        algo=unique(dat$algo),
        x=xy_vals,
        y=xy_vals
    )

    na_tab <- merge(
        full_grid,
        na_tab,
        by=c("algo", "x", "y"),
        all.x=TRUE
    )

    na_tab$x <- factor(
        na_tab$x,
        levels=xy_vals,
        labels=xy_labels
    )

    na_tab$y <- factor(
        na_tab$y,
        levels=rev(xy_vals),
        labels=rev(xy_labels)
    )

    p <- ggplot(
        na_tab,
        aes(
            x=x,
            y=y,
            fill=n_na
        )
    ) +
        geom_tile(
            colour="white",
            linewidth=0.4
        ) +
        geom_text(
            aes(
                label=ifelse(is.na(n_na), "NA", n_na)
            ),
            family=FONT_FAMILY,
            size=5.5
        ) +
        facet_wrap(
            ~algo,
            nrow=1
        ) +
        scale_fill_gradient(
            low="white",
            high="red",
            na.value="grey85",
            name="# undefined"
        ) +
        labs(
            x=expression(lambda[C]),
            y=expression(lambda[E])
        ) +
        coord_equal() +
        expanded_theme +
        theme(
            axis.text.x=element_text(
                angle=45,
                hjust=1,
                size=16
            ),
            axis.text.y=element_text(
                size=16
            ),
            axis.title=element_text(
                size=20
            ),
            strip.text=element_text(
                size=20
            ),
            legend.title=element_text(
                size=15
            ),
            legend.text=element_text(
                size=14
            )
        )

    ggsave(
        outfile,
        p,
        width=13,
        height=7,
        device=cairo_pdf
    )
}




################################################################################
# Make plots using the same grouping idea as your current final plots:
# edge retaining probability + k + expression condition
#
# Prefix/species are pooled within each group.
################################################################################

expanded_tab$group <- paste(
    expanded_tab$edge_adding_prob,
    expanded_tab$ldak,
    expanded_tab$exp_mat,
    sep="_"
)

for(g in unique(expanded_tab$group)){

    cur <- expanded_tab[expanded_tab$group == g, ]

    # Pearson
    plot_expanded_boxplot(
        cur,
        "Pearson",
        paste0(
            expanded_outdir,
            "/Pearson_boxplot_",
            g,
            ".pdf"
        )
    )

    plot_expanded_heatmap(
        cur,
        "Pearson",
        paste0(
            expanded_outdir,
            "/Pearson_heatmap_",
            g,
            ".pdf"
        )
    )

    # Spearman
    plot_expanded_boxplot(
        cur,
        "Spearman",
        paste0(
            expanded_outdir,
            "/Spearman_boxplot_",
            g,
            ".pdf"
        )
    )

    plot_expanded_heatmap(
        cur,
        "Spearman",
        paste0(
            expanded_outdir,
            "/Spearman_heatmap_",
            g,
            ".pdf"
        )
    )

    plot_expanded_NA_heatmap(
        cur,
        "Pearson",
        paste0(
            expanded_outdir,
            "/Pearson_NA_heatmap_",
            g,
            ".pdf"
        )
    )

    plot_expanded_NA_heatmap(
        cur,
        "Spearman",
        paste0(
            expanded_outdir,
            "/Spearman_NA_heatmap_",
            g,
            ".pdf"
        )
    )
}