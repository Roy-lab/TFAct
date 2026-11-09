# Input expression should be/transformed to a normalized matrix; prior can have signed.
library(TIGER)

args <- commandArgs(trailingOnly = TRUE)

in_exp=args[[1]]#"inputtoTFA/0.2/k20/Calb/rep0/E.txt"
tf_lst=args[[2]]#"/mnt/dv/wid/projects2/Roy-regnet-inference/projects/multinetlearning/multinetlearning/speciesspecsim/testinput/tflist.txt"
prior_path=args[[3]]#"inputtoTFA/0.2/k20/Calb_prior.txt"
out_A_path=args[[4]]
out_P_path=args[[5]]

raw_dat = read.table(in_exp, check.names=FALSE, sep='\t', row.names=1)
TF_used = read.table(tf_lst, check.names=FALSE)[,"V1"]
prior_tab = read.table(prior_path, check.names=FALSE, sep='\t')
colnames(raw_dat) <- paste0("cell", 1:ncol(raw_dat))

expr <- (raw_dat - apply(raw_dat, 1, mean))/apply(raw_dat, 1, sd) # needs normalized data
expr[is.na(expr)] <- 0

expr <- expr[apply(expr[,-1], 1, function(x) !all(x==0)),]


prior <- matrix(0, length(TF_used), nrow(expr))
rownames(prior) <- TF_used
colnames(prior) <- rownames(expr)

for(i in 1:nrow(prior_tab)){
    if(prior_tab[i,1] %in% rownames(prior) && prior_tab[i,2] %in% colnames(prior)){
       prior[prior_tab[i,1], prior_tab[i,2]] = 1
    }
}

ss = TIGER(expr,prior)

omitted_genes <- setdiff(rownames(raw_dat),rownames(ss$W))
omitted_TFs <- setdiff(TF_used,colnames(ss$W))

# make the new grn and tfa have the same desired dimension. If TF or genes are filtered, set the corresponding entries to 0
A <- matrix(0, nrow(raw_dat), length(TF_used))
colnames(A) <- TF_used
rownames(A) <- rownames(raw_dat)

A[rownames(ss$W), colnames(ss$W)] = ss$W

P <- matrix(0, length(TF_used), ncol(raw_dat))
rownames(P) <- TF_used
colnames(P) <- colnames(raw_dat)

P[rownames(ss$Z), colnames(ss$Z)] <- ss$Z

write.table(format(A, digits=5), out_A_path, quote = FALSE, sep = "\t",row.names = FALSE,
            col.names = FALSE)
write.table(format(P, digits=5), out_P_path, quote = FALSE, sep = "\t",row.names = FALSE,
            col.names = FALSE)