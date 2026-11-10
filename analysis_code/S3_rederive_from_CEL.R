suppressMessages({library(affy)}); source("score_lib_nofilter.R")
C <- "/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
SP <- read.csv(file.path(C,"_contrast_spec.csv"), stringsAsFactors=FALSE); sp <- SP[SP$study=="S3_GSE36878",]
an <- read.csv(file.path(C,"GPL11180_annot.tsv"), sep="\t", check.names=FALSE, quote="")
map <- setNames(as.character(an$symbol), as.character(an$ID))
fs <- list.files("/tmp/cel36878", full.names=TRUE)
gsm <- sub("_.*","",basename(fs))
pub <- list(ISG_IFNA=c(-0.000119,0.002076), ISG_IFNG=c(-0.000046,0.001771),
            ISG_ISGF3=c(-0.000223,0.002744), YAP=c(-0.000969,0.004461))
run <- function(tag, files, labs){
  ab <- suppressMessages(ReadAffy(filenames=files))
  e <- exprs(suppressMessages(rma(ab))); colnames(e) <- labs
  s <- unname(map[rownames(e)]); ok <- !is.na(s)&s!=""&s!="---"
  M <- e[ok,,drop=FALSE]; rownames(M) <- s[ok]
  o <- order(rowMeans(M), decreasing=TRUE); M <- M[o,,drop=FALSE]; X <- M[!duplicated(rownames(M)),,drop=FALSE]
  cat(sprintf("--- %s : %d probes -> %d symbols ---\n", tag, nrow(e), nrow(X)))
  for (mo in names(pub)) {
    hit <- intersect(GENESETS[[paste0(mo,"_m")]], rownames(X))
    v <- simpleScore(rankGenes(X), upSet=hit, centerScore=TRUE)$TotalScore; names(v) <- colnames(X)
    v <- v[sp$sample]; g <- ifelse(sp$group=="treated",1,0); st <- as.character(sp$stratum)
    ds <- tapply(seq_along(v), st, function(k) mean(v[k][g[k]==1]) - mean(v[k][g[k]==0]))
    d <- mean(ds); se <- sqrt(var(v[g==1])/9 + var(v[g==0])/9)
    hitflag <- if (abs(d-pub[[mo]][1])<3e-5 && abs(se-pub[[mo]][2])<3e-5) "  <<< MATCH" else ""
    cat(sprintf("  %-10s hit %3d  delta %+0.6f (pub %+0.6f)  se %0.6f (pub %0.6f)%s\n",
        mo, length(hit), d, pub[[mo]][1], se, pub[[mo]][2], hitflag)) } }
run("RMA on all 21 CELs (as the depositors did)", fs, gsm)
k <- gsm %in% sp$sample
run("RMA on the 18 contrast CELs only", fs[k], gsm[k])
