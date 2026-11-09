
library(data.table)
library(ggplot2)

###############################################################################
# Settings
###############################################################################

indir <- "plots/tradeseq"
outdir <- file.path(indir, "comparison")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

algos <- c(
    "ibpg_2fct_500_0",
    "ibpg_2fct",
    "ibpg_2fct_500_1000",
    "supirfactor_hierarchical",
    "TIGER"
)

methods <- c("DE", algos)
tests <- c("startVsEnd", "association", "diffEnd", "pattern")

FONT_FAMILY <- "Arial"

###############################################################################
# Read top genes
###############################################################################

read_top_genes <- function(method, test){

    f <- if(method == "DE"){
        file.path(indir, "DE", test, "top10_genes.txt")
    } else {
        file.path(indir, "TFA", method, test, "top10_genes.txt")
    }

    if(!file.exists(f)){
        warning("Missing: ", f)
        return(NULL)
    }

    fread(f, check.names=FALSE)
}

###############################################################################
# Pairwise Jaccard
###############################################################################

compute_jaccard <- function(test){

    tabs <- setNames(
        lapply(methods, read_top_genes, test=test),
        methods
    )

    pairs <- combn(methods, 2, simplify=FALSE)
    results <- list()

    for(pair in pairs){

        a <- pair[1]
        b <- pair[2]

        if(is.null(tabs[[a]]) || is.null(tabs[[b]]))
            next

        # Only comparisons available in both methods
        lineages <- intersect(names(tabs[[a]]), names(tabs[[b]]))

        for(lineage in lineages){

            genes_a <- unique(trimws(as.character(tabs[[a]][[lineage]])))
            genes_b <- unique(trimws(as.character(tabs[[b]][[lineage]])))

            genes_a <- genes_a[!is.na(genes_a) & nzchar(genes_a)]
            genes_b <- genes_b[!is.na(genes_b) & nzchar(genes_b)]

            if(length(genes_a) == 0 || length(genes_b) == 0)
                next

            shared <- intersect(genes_a, genes_b)
            union_genes <- union(genes_a, genes_b)

            results[[length(results)+1]] <- data.table(
                test=test,
                lineage=lineage,
                method1=a,
                method2=b,
                n1=length(genes_a),
                n2=length(genes_b),
                intersection=length(shared),
                union=length(union_genes),
                jaccard=length(shared)/length(union_genes),
                shared_genes=paste(shared, collapse=",")
            )
        }
    }

    if(length(results) == 0)
        return(NULL)

    rbindlist(results)
}

###############################################################################
# heatmap
###############################################################################


plot_jaccard <- function(dat, test){

    labels <- c(
        "DE"="DE",
        "ibpg_2fct_500_0"="TFAct-F_500_0",
        "supirfactor_hierarchical"="SupirFactor hierarchical",
        "TIGER"="TIGER",
        "ibpg_2fct"="TFAct-F_500_100",
        "ibpg_2fct_500_1000"="TFAct-F_500_1000"
    )

    method_order <- unname(labels[methods])

    # Mirror pairwise comparisons
    mat <- rbind(
        dat[, .(lineage, method1, method2, jaccard)],
        dat[, .(
            lineage,
            method1=method2,
            method2=method1,
            jaccard
        )]
    )

    # Add diagonal (only for methods with available data)
    diagonal <- unique(rbind(
        dat[, .(lineage, method=method1)],
        dat[, .(lineage, method=method2)]
    ))

    diagonal[, `:=`(
        method1=method,
        method2=method,
        jaccard=1
    )]

    mat <- rbind(
        mat,
        diagonal[, .(lineage, method1, method2, jaccard)]
    )

    # Display labels
    mat[, method1 := factor(
        labels[method1], levels=method_order
    )]

    mat[, method2 := factor(
        labels[method2], levels=rev(method_order)
    )]

    mat[, lineage := factor(
        lineage, levels=unique(dat$lineage)
    )]

    p <- ggplot(mat, aes(
        x=method1, y=method2, fill=jaccard
    )) +
        geom_tile(colour="white", linewidth=0.5) +
        geom_text(
            aes(label=sprintf("%.2f", jaccard)),
            family=FONT_FAMILY,
            size=4.2
        ) +
        facet_wrap(~lineage, ncol=3) +
        scale_fill_gradient(
            low="white",
            high="#2166AC",
            limits=c(0, 1),
            name="Jaccard"
        ) +
        coord_equal() +
        labs(x=NULL, y=NULL) +
        theme_minimal(
            base_family=FONT_FAMILY,
            base_size=18
        ) +
        theme(
            axis.text.x=element_text(
                angle=45, hjust=1,
                size=12, colour="black"
            ),
            axis.text.y=element_text(
                size=12, colour="black"
            ),
            strip.text=element_text(
                size=19, face="bold"
            ),
            panel.grid=element_blank(),
            strip.background=element_blank(),
            panel.spacing=grid::unit(1.2, "lines"),
            legend.title=element_text(size=17),
            legend.text=element_text(size=15),
            plot.margin=margin(10, 10, 10, 10)
        )

    ggsave(
        file.path(outdir, paste0(test, "_jaccard.pdf")),
        p,
        width=18,
        height=12,
        device=cairo_pdf
    )
}


###############################################################################
# Main
###############################################################################

for(test in tests){

    message("Processing: ", test)

    dat <- compute_jaccard(test)

    if(is.null(dat)){
        message("No valid comparisons for ", test)
        next
    }

    fwrite(
        dat,
        file.path(outdir, paste0(test, "_jaccard.txt")),
        sep="\t",
        quote=FALSE,
        na="NA"
    )

    plot_jaccard(dat, test)
}
