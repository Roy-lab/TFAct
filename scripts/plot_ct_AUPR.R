library(ggplot2)
library(parallel)

###############################################################################
# Settings
###############################################################################

FONT_FAMILY <- "Arial"

outdir <- "TFA_output/analysis/ct/AUPR_expanded"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

algos <- c("ibpg_2fct", "kl_nmf_2fct")

Prefixes <- c(
    "Spom", "Anc2", "Anc5", "Calb",
    "Cgla", "right", "Sbay"
)

Exp_mat_w_dropouts <- c("E", "E_0.6")

Edge_adding_probs_expanded <- c(0.4)
Ldaks_expanded <- c(20)
TFA_sim_reps <- c("rep0")

xy_vals <- c(
    "0", "10", "100", "1000",
    "10000", "100000", "1000000", "10000000"
)

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

# Four combinations:
# AUPR/AUPR_TF.txt
# AUPR/AUPR.txt
# AUPR_prior/AUPR_TF.txt
# AUPR_prior/AUPR.txt
plot_sets <- expand.grid(
    folder=c("AUPR", "AUPR_prior"),
    file=c("AUPR_TF.txt", "AUPR.txt"),
    stringsAsFactors=FALSE
)


###############################################################################
# Theme
###############################################################################

expanded_theme <- theme_classic(
    base_family=FONT_FAMILY,
    base_size=18
) +
    theme(
        axis.text=element_text(size=16, colour="black"),
        axis.title=element_text(size=20),
        strip.text=element_text(size=18),
        strip.background=element_blank(),
        panel.spacing.x=grid::unit(0.3, "lines"),
        plot.margin=margin(8, 8, 8, 8)
    )


###############################################################################
# Read one AUPR file
###############################################################################

read_AUPR_values <- function(f){

    dat <- read.table(
        f,
        header=FALSE,
        sep="",
        stringsAsFactors=FALSE,
        fill=TRUE
    )

    # Find numeric columns
    numeric_versions <- lapply(
        dat,
        function(z) suppressWarnings(as.numeric(z))
    )

    n_numeric <- sapply(
        numeric_versions,
        function(z) sum(is.finite(z))
    )

    if(max(n_numeric) == 0){
        return(numeric(0))
    }

    # use column with most numeric values
    vals <- numeric_versions[[which.max(n_numeric)]]

    return(vals)
}


###############################################################################
# Collect expanded AUPR data
###############################################################################

read_expanded_AUPR <- function(
    edge_adding_prob,
    ldak,
    prefix,
    exp_mat,
    sim_rep,
    aupr_folder,
    aupr_file
){

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
                    "_stability/",
                    aupr_folder, "/",
                    aupr_file
                )

                if(!file.exists(cur_f)){
                    next
                }

                vals <- read_AUPR_values(cur_f)

                if(length(vals) == 0){
                    next
                }

                cur_df <- data.frame(
                    value=vals,
                    algo=ifelse(
                        algo == "ibpg_2fct",
                        "F",
                        "KL"
                    ),
                    x=x,
                    y=y,
                    prefix=prefix,
                    stringsAsFactors=FALSE
                )

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


###############################################################################
# Boxplots
###############################################################################

plot_expanded_AUPR_boxplot <- function(dat, outfile){

    dat$x <- factor(
        dat$x,
        levels=xy_vals
    )

    dat$y <- factor(
        dat$y,
        levels=xy_vals
    )

    p <- ggplot(
        dat,
        aes(
            x=y,
            y=value
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
                x=as_labeller(
                    function(z)
                        paste0(
                            "lambda[C] == '",
                            xy_labels[z],
                            "'"
                        ),
                    label_parsed
                )
            )
        ) +
        scale_x_discrete(
            labels=xy_labels,
            drop=FALSE
        ) +
        labs(
            x=expression(lambda[E]),
            y="AUPR"
        ) +
        expanded_theme +
        theme(
            axis.text.x=element_text(
                angle=45,
                hjust=1,
                vjust=1,
                size=11
            ),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=20),
            strip.text=element_text(size=17)
        )

    ggsave(
        outfile,
        p,
        width=14,
        height=7,
        device=cairo_pdf
    )
}


###############################################################################
# Median heatmaps
###############################################################################

plot_expanded_AUPR_heatmap <- function(dat, outfile){

    med_tab <- aggregate(
        dat$value,
        by=list(
            algo=dat$algo,
            x=dat$x,
            y=dat$y
        ),
        FUN=function(z){

            z <- z[is.finite(z)]

            if(length(z) == 0){
                NA_real_
            } else {
                median(z)
            }
        }
    )

    colnames(med_tab)[4] <- "value"

    full_grid <- expand.grid(
        algo=unique(dat$algo),
        x=xy_vals,
        y=xy_vals,
        stringsAsFactors=FALSE
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
            linewidth=0.4
        ) +
        geom_text(
            aes(
                label=ifelse(
                    is.na(value),
                    "NA",
                    sprintf("%.2f", value)
                )
            ),
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
            midpoint=0.5,
            limits=c(0, 1),
            na.value="grey90",
            name="Median\nAUPR"
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
                size=15
            ),
            axis.text.y=element_text(size=16),
            axis.title=element_text(size=20),
            strip.text=element_text(size=20),
            legend.title=element_text(size=15),
            legend.text=element_text(size=14)
        )

    ggsave(
        outfile,
        p,
        width=13,
        height=7,
        device=cairo_pdf
    )
}


###############################################################################
# Main
###############################################################################

for(exp_mat in Exp_mat_w_dropouts){

    for(k in seq_len(nrow(plot_sets))){

        aupr_folder <- plot_sets$folder[k]
        aupr_file <- plot_sets$file[k]

        message(
            "Processing ",
            aupr_folder,
            "/",
            aupr_file,
            "  ",
            exp_mat
        )

        stat_lst <- mclapply(
            Prefixes,
            function(prefix){

                read_expanded_AUPR(
                    edge_adding_prob=0.4,
                    ldak=20,
                    prefix=prefix,
                    exp_mat=exp_mat,
                    sim_rep="rep0",
                    aupr_folder=aupr_folder,
                    aupr_file=aupr_file
                )
            },
            mc.cores=7
        )

        stat_lst <- stat_lst[
            !sapply(stat_lst, is.null)
        ]

        if(length(stat_lst) == 0){
            message("No files found; skipping.")
            next
        }

        dat <- do.call(
            rbind,
            stat_lst
        )

        #######################################################################
        # Filename without redundant AUPR_ prefix
        #######################################################################

        file_tag <- sub(
            "\\.txt$",
            "",
            aupr_file
        )

        if(aupr_folder == "AUPR_prior"){
            file_tag <- paste0(
                "prior_",
                file_tag
            )
        }

        suffix <- paste0(
            file_tag,
            "_",
            exp_mat
        )

        plot_expanded_AUPR_boxplot(
            dat,
            paste0(
                outdir,
                "/",
                suffix,
                "_boxplot.pdf"
            )
        )

        plot_expanded_AUPR_heatmap(
            dat,
            paste0(
                outdir,
                "/",
                suffix,
                "_heatmap.pdf"
            )
        )
    }
}