
source("score_lib.R")
suppressMessages({library(org.Mm.eg.db); library(org.Hs.eg.db); library(readxl)})
dir.create("out", showWarnings=FALSE)

ens2sym <- function(ids, sp) {
  ids <- sub("\\..*$", "", ids)
  db <- if (sp=="mouse") org.Mm.eg.db else org.Hs.eg.db
  unname(suppressMessages(AnnotationDbi::mapIds(db, keys=ids, keytype="ENSEMBL", column="SYMBOL", multiVals="first")))
}
by_sym <- function(mat, sym) { k <- !is.na(sym) & sym!=""; mat <- mat[k,,drop=FALSE]; rownames(mat) <- sym[k]; mat }
rd <- function(p, ...) read.delim(p, check.names=FALSE, ...)

L <- list()
L$GSE304867 <- function() { m <- read.csv("data/GSE304867/GSE304867_count_matrix.csv.gz", row.names=1, check.names=FALSE)
  m <- by_sym(as.matrix(m), ens2sym(rownames(m),"mouse"))
  list(expr=prep_expr(m,"counts"), grp=c(rep(0,6),rep(1,6)), sp="mouse") }
L$GSE285024 <- function() { d <- rd("data/GSE285024/GSE285024_Normalized_Count_Matrix_BMDM_confiner.txt.gz", sep=" ")
  cn <- c("WT_NC1","WT_NC3","WT_3_2","WT_3_3","WT_3_6"); m <- as.matrix(d[,cn]); rownames(m) <- d$Symbol
  list(expr=prep_expr(m,"norm"), grp=c(0,0,1,1,1), sp="mouse") }
L$GSE303278 <- function() { m <- as.matrix(rd("data/GSE303278/GSE303278_on_vs_off_results_FPKM_DESeq.tsv.gz", row.names=1))
  list(expr=prep_expr(m,"fpkm"), grp=c(0,0,0,1,1,1), sp="mouse") }
L$GSE344008 <- function() { d <- rd("data/GSE344008/GSE344008_gene_expression.txt.gz")
  cc <- grep("^Count ", names(d), value=TRUE); m <- as.matrix(d[,cc]); rownames(m) <- d[["Gene Name"]]
  list(expr=prep_expr(m,"counts"), grp=c(0,0,0,1,1,1), sp="mouse") }
L$GSE262146 <- function() { fs <- sort(list.files("data/GSE262146", pattern="count.txt.gz$", full.names=TRUE))
  ms <- lapply(fs, function(f){ d <- rd(f, skip=1); setNames(d[[7]], d$Geneid) })
  g <- Reduce(intersect, lapply(ms, names)); m <- sapply(ms, function(x) x[g]); rownames(m) <- g
  m <- by_sym(m, ens2sym(rownames(m),"human"))
  list(expr=prep_expr(m,"counts"), grp=c(0,0,0,1,1,1), sp="human") }
L[["E-MTAB-15203"]] <- function() { d <- rd("data/E-MTAB-15203_raw_counts_table.tsv", row.names=1)
  m <- as.matrix(d); m <- by_sym(m, ens2sym(rownames(m),"human"))
  m <- m[, c(paste0("D",1:6,"_DMSO"), paste0("D",1:6,"_Y1"))]
  list(expr=prep_expr(m,"counts"), grp=c(rep(0,6),rep(1,6)), sp="human", paired=rep(paste0("D",1:6),2)) }
L$GSE290857 <- function() { d <- rd("data/GSE290857/GSE290857_1_genes_fpkm_expression.txt.gz")
  cc <- grep("^FPKM[.]", names(d), value=TRUE); m <- as.matrix(d[,cc]); rownames(m) <- d$gene_name
  list(expr=prep_expr(m,"fpkm"), grp=c(rep(0,4),rep(1,4)), sp="mouse") }
L$GSE331328 <- function() { d <- read.csv("data/GSE331328/GSE331328_gene_sample_FPKM.csv.gz", row.names=1, check.names=FALSE)
  m <- by_sym(as.matrix(d), ens2sym(rownames(d),"mouse"))
  list(expr=prep_expr(m,"fpkm"), grp=c(0,0,0,1,1,1), sp="mouse") }
L$GSE273466 <- function() { d <- rd("data/GSE273466/GSE273466_gene_exp.xls.gz")
  cc <- c("LPS-1_COUNT","LPS-2_COUNT","LPS-3_COUNT","LIPUS-LPS-1_COUNT","LIPUS-LPS-2_COUNT","LIPUS-LPS-3_COUNT")
  m <- as.matrix(d[,cc]); rownames(m) <- d$ID
  list(expr=prep_expr(m,"counts"), grp=c(0,0,0,1,1,1), sp="mouse") }
L$GSE314256 <- function() { fs <- sort(list.files("data/GSE314256", pattern="quant[.]genes[.]sf[.]gz$", full.names=TRUE))
  ms <- lapply(fs, function(f){ d <- rd(f); tapply(d$TPM, d$GeneName, max) })
  g <- Reduce(intersect, lapply(ms, names)); m <- sapply(ms, function(x) x[g]); rownames(m) <- g
  list(expr=prep_expr(m,"tpm"), grp=c(0,0,0,1,1,1), sp="human") }
L$GSE236912 <- function() { d <- rd("data/GSE236912/GSE236912_all.fpkm_anno.xls.gz")
  cc <- c(paste0("m1-",1:4), paste0("m3-",1:4)); m <- as.matrix(d[,cc]); rownames(m) <- d$GeneSymbol
  list(expr=prep_expr(m,"counts"), grp=c(rep(0,4),rep(1,4)), sp="mouse") }
L$GSE273977 <- function() { fs <- sort(list.files("data/GSE273977", pattern="readcount-sample", full.names=TRUE))
  ms <- lapply(fs, function(f){ d <- rd(f); setNames(d[[2]], d$GeneID) })
  g <- Reduce(intersect, lapply(ms, names)); m <- sapply(ms, function(x) x[g]); rownames(m) <- g
  list(expr=prep_expr(m,"counts"), grp=c(0,0,0,1,1,1), sp="mouse") }
L$GSE133069 <- function() { d <- as.data.frame(read_excel("data/GSE133069/GSE133069_Gene_expression_WT_KO_Cyclic_Static.xlsx"))
  m <- as.matrix(d[,c("WT_Static_NT_8","WT_Cyclic_NT_15")]); rownames(m) <- d$Gene_Name
  list(expr=prep_expr(m,"fpkm"), grp=c(0,1), sp="mouse") }
L$GSE190470 <- function() { d <- rd("data/GSE190470/GSE190470_gene_count_matrix.txt.gz", row.names=1)
  m <- by_sym(as.matrix(d), ens2sym(rownames(d),"mouse"))
  list(expr=prep_expr(m,"counts"), grp=c(0,0,0,1,1,1), sp="mouse") }
L$GSE125806 <- function() { d <- rd("data/GSE125806/GSE125806_RAWCOUNT_MATRIX.txt.gz", row.names=1)
  list(expr=prep_expr(as.matrix(d),"counts"), grp=c(0,1), sp="mouse") }
L$GSE185761 <- function() { d <- rd("data/GSE185761/GSE185761_genecounts.txt.gz", row.names=1)
  cc <- c("Macro_Ctrl_repA","Macro_Ctrl_repB","Macro_Day14_repA","Macro_Day14_repB","Macro_Day14_repC")
  m <- by_sym(as.matrix(d[,cc]), ens2sym(rownames(d),"mouse"))
  list(expr=prep_expr(m,"counts"), grp=c(0,0,1,1,1), sp="mouse") }
L$GSE143780 <- function() { d <- rd("data/GSE143780/GSE143780_norm_gene_counts.txt.gz")
  cn <- names(d); cpb <- grep("CPB[.]00h|CPB[.]end", cn, value=TRUE)
  pid <- sub("[.]CPB.*","",cpb); tp <- ifelse(grepl("00h",cpb),0,1)
  o <- order(tp, pid); cpb <- cpb[o]
  m <- as.matrix(d[,cpb]); rownames(m) <- d$name
  list(expr=prep_expr(m,"norm"), grp=tp[o], sp="human", paired=pid[o]) }

META <- read.csv("meta_bulk.csv", stringsAsFactors=FALSE)
res <- list()
for (i in seq_len(nrow(META))) {
  id <- META$id[i]
  r <- try({ x <- L[[id]]()
    rr <- score_dataset(x$expr, x$grp, x$sp, strata=x$strata, paired_id=x$paired, nrand=1000, seed=42)
    cbind(META[i,], rr, n_ctrl=sum(x$grp==0), n_treat=sum(x$grp==1), n_genes_bg=nrow(x$expr)) }, silent=TRUE)
  if (inherits(r,"try-error")) { cat("FAIL", id, as.character(r), "\n"); next }
  res[[id]] <- r
  cat("ok", id, "ISGF3 d=", round(r$delta[r$module=="ISG_ISGF3"],4), "p_perm=", r$p_perm[r$module=="ISG_ISGF3"],
      "z=", round(r$z_rand[r$module=="ISG_ISGF3"],2), "\n")
  write.csv(do.call(rbind, res), "out/bulk_scores.csv", row.names=FALSE)
}
cat("BULKDONE\n")

