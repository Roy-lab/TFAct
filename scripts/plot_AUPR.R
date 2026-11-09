
library(ggplot2)
library(parallel)

###############################################################################
# Settings
###############################################################################

FONT_FAMILY <- "Arial"

outdir <- "TFA_output/analysis/aupr_comparison"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

algos <- c(
    "baselines/hals_prior", "ibpg", "ibpg_wogr", "ibpg_noreg",
    "nmqf_hals", "nmqf_hals_noreg", "FastNCA", "NINCA_quad",
    "pmf_grn", "regularized_nca", "BITFAM", "TIGER",
    "NetREX_CF", "supirfactor_shallow", "supirfactor_hierarchical",
    "ibpg_2f", "ibpg_2fct", "kl_nmf_2fct"
)

ignored_algos <- c(
    "ibpg_2f", "ibpg_noreg", "ibpg_wogr",
    "nmqf_hals", "nmqf_hals_noreg", "ibpg", "hals_prior"
)

algos <- algos[!basename(algos) %in% ignored_algos]

Edge_adding_probs <- c(0.2, 0.4)
Ldaks <- c(10, 20)
Prefixes <- c("Spom", "Anc2", "Anc5", "Calb", "Cgla", "right", "Sbay")
Exp_mat_w_dropouts <- c("E", "E_0.6")
TFA_sim_reps <- c("rep0")

plot_sets <- expand.grid(
    folder=c("AUPR", "AUPR_prior"),
    file=c("AUPR.txt", "AUPR_TF.txt"),
    stringsAsFactors=FALSE
)

###############################################################################
# Algorithm names (preserve original order)
###############################################################################

formal_algo_names <- function(x){

    labels <- c(
        "ibpg"="IBPG-p",
        "hals_prior"="HALS-p",
        "ibpg_2fct"="TFAct-F_500_100",
        "kl_nmf_2fct"="TFAct-KL_500_100"
    )

    x <- basename(x)
    idx <- x %in% names(labels)
    x[idx] <- labels[x[idx]]
    x
}

algo_order <- unique(formal_algo_names(algos))

###############################################################################
# Read AUPR values
###############################################################################

read_AUPR_values <- function(f){

    dat <- tryCatch(
        read.table(
            f, header=FALSE, sep="",
            stringsAsFactors=FALSE, fill=TRUE
        ),
        error=function(e) NULL
    )

    if(is.null(dat) || ncol(dat) == 0)
        return(numeric(0))

    numeric_versions <- lapply(
        dat, function(z) suppressWarnings(as.numeric(z))
    )

    n_numeric <- sapply(
        numeric_versions, function(z) sum(is.finite(z))
    )

    if(max(n_numeric) == 0)
        return(numeric(0))

    vals <- numeric_versions[[which.max(n_numeric)]]
    vals[is.finite(vals)]
}

###############################################################################
# Collect all conditions and algorithms
###############################################################################

read_AUPR_all <- function(aupr_folder, aupr_file){

    settings <- expand.grid(
        algo=algos,
        edge=Edge_adding_probs,
        k=Ldaks,
        exp_mat=Exp_mat_w_dropouts,
        prefix=Prefixes,
        sim_rep=TFA_sim_reps,
        stringsAsFactors=FALSE
    )

    dfs <- mclapply(seq_len(nrow(settings)), function(i){

        s <- settings[i, ]

        cur_f <- paste0(
            "TFA_output/", s$algo, "/",
            s$edge, "/k", s$k, "/",
            s$prefix, "/", s$sim_rep, "/",
            s$exp_mat, "_stability/",
            aupr_folder, "/", aupr_file
        )

        if(!file.exists(cur_f))
            return(NULL)

        vals <- read_AUPR_values(cur_f)

        if(length(vals) == 0)
            return(NULL)

        data.frame(
            algo=formal_algo_names(s$algo),
            prefix=s$prefix,
            edge=s$edge,
            k=s$k,
            exp_mat=s$exp_mat,
            sim_rep=s$sim_rep,
            value=vals,
            stringsAsFactors=FALSE
        )
    }, mc.cores=8)

    dfs <- Filter(Negate(is.null), dfs)

    if(length(dfs) == 0)
        return(NULL)

    dat <- do.call(rbind, dfs)

    dat$algo <- factor(dat$algo, levels=algo_order)

    dat$density <- factor(
        ifelse(dat$edge == 0.4, "Denser", "Sparser"),
        levels=c("Sparser", "Denser")
    )

    dat$noise <- factor(
        ifelse(dat$k == 20, "Lower noise", "Higher noise"),
        levels=c("Higher noise", "Lower noise")
    )

    dat$dropout <- factor(
        ifelse(dat$exp_mat == "E", "No dropout", "Dropout 0.6"),
        levels=c("No dropout", "Dropout 0.6")
    )

    dat
}

###############################################################################
# Shared theme
###############################################################################

aupr_theme <- theme_classic(
    base_family=FONT_FAMILY,
    base_size=18
) +
    theme(
        axis.text.x=element_text(
            angle=45, hjust=1, vjust=1,
            colour="black"
        ),
        axis.text.y=element_text(colour="black"),
        axis.title.y=element_text(margin=margin(r=10)),
        axis.line=element_line(linewidth=0.6),
        axis.ticks=element_line(linewidth=0.5),
        strip.background=element_blank(),
        plot.margin=margin(10, 12, 10, 10)
    )

###############################################################################
# Individual boxplots
###############################################################################

plot_AUPR_individual <- function(dat, outfile){

    p <- ggplot(dat, aes(x=algo, y=value)) +
        geom_boxplot(
            width=0.65,
            fill="#A9C9E2",
            colour="#253746",
            linewidth=0.55,
            outlier.size=0.8,
            outlier.alpha=0.45,
            na.rm=TRUE
        ) +
        scale_x_discrete(drop=TRUE) +
        coord_cartesian(ylim=c(0, 1)) +
        labs(x=NULL, y="AUPR") +
        aupr_theme +
        theme(
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=17),
            axis.title.y=element_text(size=21)
        )

    ggsave(
        outfile, p,
        width=12, height=6.5,
        device=cairo_pdf
    )
}

###############################################################################
# Combined boxplots: 4 rows x 2 columns
###############################################################################

plot_AUPR_combined <- function(dat, outfile){

    p <- ggplot(dat, aes(x=algo, y=value)) +
        geom_boxplot(
            width=0.65,
            fill="#A9C9E2",
            colour="#253746",
            linewidth=0.5,
            outlier.size=0.6,
            outlier.alpha=0.4,
            na.rm=TRUE
        ) +
        facet_grid(
            rows=vars(noise, density),
            cols=vars(dropout),
            switch="y"
        ) +
        scale_x_discrete(drop=TRUE) +
        coord_cartesian(ylim=c(0, 1)) +
        labs(x=NULL, y="AUPR") +
        aupr_theme +
        theme(
            axis.text.x=element_text(size=13, angle=45, hjust=1),
            axis.text.y=element_text(size=17),
            axis.title.y=element_text(size=21),
            strip.text=element_text(size=18, colour="black"),
            panel.spacing=grid::unit(0.4, "lines"),
            strip.placement="outside"
        )

    ggsave(
        outfile, p,
        width=9, height=12,
        device=cairo_pdf
    )
}

###############################################################################
# Main: individual + combined
###############################################################################

settings <- expand.grid(
    edge=Edge_adding_probs,
    k=Ldaks,
    exp_mat=Exp_mat_w_dropouts,
    sim_rep=TFA_sim_reps,
    stringsAsFactors=FALSE
)

for(j in seq_len(nrow(plot_sets))){

    aupr_folder <- plot_sets$folder[j]
    aupr_file <- plot_sets$file[j]

    file_tag <- sub("\\.txt$", "", aupr_file)

    if(aupr_folder == "AUPR_prior")
        file_tag <- paste0("prior_", file_tag)

    message("Processing: ", file_tag)

    # Read once and reuse for every plot
    dat <- read_AUPR_all(aupr_folder, aupr_file)

    if(is.null(dat)){
        message("No data found; skipping.")
        next
    }

    ###########################################################################
    # Individual plots
    ###########################################################################

    for(i in seq_len(nrow(settings))){

        s <- settings[i, ]

        subset_dat <- dat[
            dat$edge == s$edge &
            dat$k == s$k &
            dat$exp_mat == s$exp_mat &
            dat$sim_rep == s$sim_rep,
        ]

        if(nrow(subset_dat) == 0)
            next

        suffix <- paste0(
            s$edge, "_k", s$k, "_",
            s$exp_mat, "_", s$sim_rep
        )

        outfile <- file.path(
            outdir,
            paste0(file_tag, "_", suffix, ".pdf")
        )

        plot_AUPR_individual(subset_dat, outfile)
    }

    ###########################################################################
    # Combined plot
    ###########################################################################

    outfile <- file.path(
        outdir,
        paste0(file_tag, "_combined.pdf")
    )

    plot_AUPR_combined(dat, outfile)

    message("Saved combined: ", outfile)
}
