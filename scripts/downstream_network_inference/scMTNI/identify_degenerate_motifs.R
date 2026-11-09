library(parallel)
library(igraph)
library(ggplot2)
library(data.table)
args=commandArgs(trailingOnly = T)
indir <- args[1] # where prior networks are located "../real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/"
cellfile <- args[2]
suffix <- args[3] # "_network.txt"
cutoff <- args[4] # 0.995
datadir <- args[5] # where cluster_genes_filteredlowexpression.txt are stored. ../real_dat_validation/scMTNI_data/scMTNI_sourcedata/processed_FBS_reprogramming/inputdata/scMTNI/data/ibpg_2fct/
outdir <- args[6] # 

theme_set(
    theme_classic(base_size = 16, base_family = "Arial")
)

cells=read.table(cellfile,stringsAsFactors = F)
cells=cells$V1

files <- file.path(indir, paste0(cells, suffix))

allres <- do.call(rbind, mclapply(files, function(f) {
    cl <- sub(paste0(suffix, "$"), "", basename(f))
    d <- read.table(f, sep = "\t", header = FALSE,
                    stringsAsFactors = FALSE)
    colnames(d) <- c("TF", "target", "score")
    tfs <- unique(d$TF)
    targets <- unique(d$target)
    M <- matrix(
        0,
        nrow = length(tfs),
        ncol = length(targets),
        dimnames = list(tfs, targets)
    )
    M[cbind(
        match(d$TF, tfs),
        match(d$target, targets)
    )] <- d$score

    ## pairwise Pearson correlations
    R <- suppressWarnings(cor(t(M), method = "pearson"))
    ## shared / union target counts
    B <- M > 0
    Nshared <- B %*% t(B)
    Ntarget <- rowSums(B)
    Nunion <- outer(Ntarget, Ntarget, "+") - Nshared

    ij <- which(upper.tri(R), arr.ind = TRUE)

    data.frame(
        cluster = cl,
        TF1 = rownames(R)[ij[, 1]],
        TF2 = colnames(R)[ij[, 2]],
        n_union = Nunion[ij],
        n_shared = Nshared[ij],
        pearson = R[ij]
    )

}, mc.cores = 20))

allres <- allres[order(-allres$pearson), ]

allres$cluster_plot <- sub("^cluster([0-9]+)$", "C\\1", allres$cluster)

# histogram for motif degeneracy
p <- ggplot(allres, aes(x = pearson)) +
    geom_histogram(
        bins = 50,
        color = "black",
        fill = "grey75",
        linewidth = 0.3
    ) +
    scale_y_log10() +
    facet_wrap(~ cluster_plot) +
    labs(
        x = "Pearson correlation",
        y = "Number of TF pairs (log10 scale)"
    ) +
    theme_classic(base_size = 16) +
    theme(
        strip.text = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 14)
    )

ggsave(
    file.path(outdir, "motif_pairwise_pearson_histograms.pdf"),
    p,
    width = 10,
    height = 7,
    device = cairo_pdf
)

# find connected components in the TF-by-TF graph where edges whose scores are >= cutoff are retained
get_families <- function(x, cutoff = 0.995) {
    x <- x[x$pearson >= cutoff, c("TF1", "TF2")]
    g <- graph_from_data_frame(x, directed = FALSE)
    comp <- components(g)$membership
    fam <- split(names(comp), comp)
    data.frame(
        family = paste0("family", seq_along(fam)),
        TFs = sapply(fam, paste, collapse = ", "),
        n_TFs = lengths(fam),
        row.names = NULL
    )
}

if (datadir != "NA"){
for (cell in cells) {
    ## motif-degenerate families
    fam <- get_families(allres[allres$cluster == cell, ])

    fwrite(fam, file.path(outdir, paste0(cell, "_motif_families.txt")), sep = "\t")

    geneskp <- fread(file.path(datadir, paste0(cell, "_genes_filteredlowexpression.txt")), header = FALSE)$V1
    geneskp <- geneskp[!grepl("_TFA$", geneskp)]

    ## retain only expressed TFs within each family
    fam_filtered <- do.call(rbind, lapply(seq_len(nrow(fam)), function(i) {
        tf <- trimws(unlist(strsplit(fam$TFs[i], ",")))
        tf <- tf[tf %in% geneskp]

        if (length(tf) <= 1)
            return(NULL)

        data.frame(family = fam$family[i], TFs = paste(tf, collapse = ", "), n_TFs = length(tf))
    }))
    fwrite(
        fam_filtered,
        file.path(outdir, paste0(cell, "_motif_families_filteredexpression.txt")),
        sep = "\t"
    )
}
}


# clusters <- cells
# nclust <- length(clusters)

# ## keep degenerate pairs
# x <- allres[allres$pearson >= cutoff, c("cluster", "TF1", "TF2")]

# ## make pair orientation consistent
# x$a <- pmin(x$TF1, x$TF2)
# x$b <- pmax(x$TF1, x$TF2)

# ## count in how many clusters each pair occurs
# pair_count <- aggregate(
#     cluster ~ a + b,
#     data = unique(x[, c("cluster", "a", "b")]),
#     FUN = length
# )
# colnames(pair_count)[3] <- "n_clusters"

# pairs_to_families <- function(x) {
#     g <- graph_from_data_frame(x[, c("a", "b")], directed = FALSE)
#     fam <- split(names(components(g)$membership),
#                  components(g)$membership)

#     data.frame(
#         family = paste0("family", seq_along(fam)),
#         TFs = sapply(fam, function(z) paste(sort(z), collapse = ", ")),
#         n_TFs = lengths(fam),
#         row.names = NULL
#     )
# }

# ## UNION: pair is degenerate in >=1 cluster
# fam_union <- pairs_to_families(pair_count)

# ## INTERSECTION: pair is degenerate in every cluster
# fam_intersection <- pairs_to_families(
#     pair_count[pair_count$n_clusters == nclust, ]
# )

# union_TFs <- sort(unique(
#     trimws(unlist(strsplit(fam_union$TFs, ",")))
# ))

# intersection_TFs <- sort(unique(
#     trimws(unlist(strsplit(fam_intersection$TFs, ",")))
# ))
