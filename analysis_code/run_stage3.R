
source("score_lib.R")
suppressMessages({library(org.Hs.eg.db); library(GEOquery)})
rd <- function(p, ...) read.delim(p, check.names=FALSE, ...)
by_sym <- function(mat, sym) { k <- !is.na(sym) & sym!=""; mat <- mat[k,,drop=FALSE]; rownames(mat) <- sym[k]; mat }

ti <- function(file, treat_keep, ti_contrast) {
  g <- getGEO(filename=file, getGPL=FALSE); e <- Biobase::exprs(g); pd <- Biobase::pData(g)
  don <- sub(".*donornr: ","",pd$characteristics_ch1.1); trt <- sub(".*treatment: ","",pd$characteristics_ch1.2)
  tiv <- sub(".*titanium: ","",pd$characteristics_ch1.3)
  s <- unname(suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys=sub("_at$","",rownames(e)),
       keytype="ENTREZID", column="SYMBOL", multiVals="first")))
  e <- by_sym(e, s)
  if (ti_contrast) { k <- trt==treat_keep; grp <- ifelse(tiv[k]=="Ti",1,0) } else {
    k <- tiv=="ctrl" & trt %in% c("NO", treat_keep); grp <- ifelse(trt[k]==treat_keep,1,0) }
  o <- order(grp, don[k]); kk <- which(k)[o]
  list(expr=prep_expr(e[,kk],"log2"), grp=grp[o], sp="human", paired=don[kk]) }

L <- list()
L$GSE99445_Ti   <- function() ti("data/pheno/GSE99445_series_matrix.txt.gz","NO",TRUE)
L$GSE99445_IFNg <- function() ti("data/pheno/GSE99445_series_matrix.txt.gz","IFN",FALSE)
L$GSE99445_IL4  <- function() ti("data/pheno/GSE99445_series_matrix.txt.gz","IL4",FALSE)
L$GSE143780 <- function() { d <- rd("data/GSE143780/GSE143780_norm_gene_counts.txt.gz")
  cpb <- grep("CPB[.]00h|CPB[.]end", names(d), value=TRUE)
  pid <- sub("[.]CPB.*","",cpb); tp <- ifelse(grepl("00h",cpb),0,1); o <- order(tp,pid)
  m <- as.matrix(d[,cpb[o]]); rownames(m) <- d$name
  list(expr=prep_expr(m,"norm"), grp=tp[o], sp="human", paired=pid[o]) }

META <- read.csv("meta_stage3.csv", stringsAsFactors=FALSE)
res <- list()
for (i in seq_len(nrow(META))) { id <- META$id[i]
  r <- try({ x <- L[[id]]()
    cat(id, "n=", length(x$grp), "pairs=", length(unique(x$paired)), "\n")
    rr <- score_dataset(x$expr, x$grp, x$sp, paired_id=x$paired, nrand=1000, seed=42)
    cbind(META[i,], rr, n_ctrl=sum(x$grp==0), n_treat=sum(x$grp==1), n_genes_bg=nrow(x$expr)) }, silent=TRUE)
  if (inherits(r,"try-error")) { cat("FAIL", id, as.character(r), "\n"); next }
  res[[id]] <- r
  cat("ok", id, "ISGF3 d=", round(r$delta[r$module=="ISG_ISGF3"],4), "p=", r$p_perm[r$module=="ISG_ISGF3"],
      "z=", round(r$z_rand[r$module=="ISG_ISGF3"],2), "\n")
  write.csv(do.call(rbind,res), "out/stage3_scores.csv", row.names=FALSE) }
cat("STAGE3DONE\n")

