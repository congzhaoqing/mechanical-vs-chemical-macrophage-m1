source("score_lib_nofilter.R")
rdm <- function(p){ d <- read.csv(p, check.names=FALSE)
  nm <- as.character(d[[1]]); m <- as.matrix(d[,-1,drop=FALSE]); rownames(m) <- nm
  m }   # duplicate symbols kept; prep_expr collapses by highest row mean
P <- rdm("out/GSE275018_pseudobulk_myeloid.csv")
cat("GSE275018 cols:", colnames(P), "\n")
arm <- substr(colnames(P),1,1)
C <- list(
 GSE275018_MNM_PBS = list(m=P[, arm %in% c("A","B")], g=ifelse(arm[arm %in% c("A","B")]=="A",0,1), s=NULL),
 GSE275018_RMF_MNM = list(m=P[, arm %in% c("B","C")], g=ifelse(arm[arm %in% c("B","C")]=="B",0,1), s=NULL),
 GSE275018_RMF_PBS = list(m=P[, arm %in% c("A","C")], g=ifelse(arm[arm %in% c("A","C")]=="A",0,1), s=NULL))
V <- rdm("out/GSE274828_pseudobulk_myeloid.csv")
V <- V[, c("ANV","AV","YNV","YV")]
C$GSE274828 <- list(m=V, g=c(0,1,0,1), s=c("aged","aged","young","young"))
B <- rdm("out/GSE186185_pseudobulk_all.csv")
B <- B[, c("sample2","sample3")]          # sample2 = Control, sample3 = Case (force)
C$GSE186185 <- list(m=B, g=c(0,1), s=NULL)
res <- list()
for (id in names(C)) {
  x <- C[[id]]
  e <- prep_expr(x$m, "counts")
  rr <- score_dataset(e, x$g, "mouse", strata=x$s, paired_id=NULL, nrand=1000, seed=42)
  res[[id]] <- cbind(id=id, rr, n_ctrl=sum(x$g==0), n_treat=sum(x$g==1), n_genes_bg=nrow(e))
  cat("ok", id, "bg=", nrow(e), "ISGF3 d=", round(rr$delta[rr$module=="ISG_ISGF3"],5),
      "TNFA d=", round(rr$delta[rr$module=="EFF_TNFA_NFKB"],5), "\n")
}
write.csv(do.call(rbind,res), "out/sc_scores_nofilter.csv", row.names=FALSE)
cat("SCSCOREDONE\n")
