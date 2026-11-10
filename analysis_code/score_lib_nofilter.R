
suppressMessages({library(singscore); library(limma); library(GSEABase)})

GSFILE <- "/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices/genesets.csv"
.gs <- read.csv(GSFILE, stringsAsFactors=FALSE)
GENESETS <- setNames(lapply(.gs$genes, function(s) unique(strsplit(s, ";")[[1]])), .gs$set)
.ef <- read.csv(file.path("/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices","genesets_effector.csv"), stringsAsFactors=FALSE)
names(.ef)[1] <- "set"
GENESETS <- c(GENESETS, setNames(lapply(.ef$genes, function(s) unique(strsplit(s, ";")[[1]])), .ef$set))
MODULES <- c("ISG_IFNA","ISG_IFNG","ISG_ISGF3","YAP","EFF_TNFA_NFKB","EFF_INFLAMM","EFF_FIELD_PANEL")

get_sets <- function(species) {
  sfx <- if (species == "human") "_h" else "_m"
  setNames(lapply(MODULES, function(m) GENESETS[[paste0(m, sfx)]]), MODULES)
}

# --- expression prep: symbols x samples, log2, duplicates collapsed, low expr dropped
prep_expr <- function(mat, scale_type) {
  mat <- as.matrix(mat); mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  if (scale_type == "counts") {
    cs <- colSums(mat); cs[cs == 0] <- 1
    mat <- log2(t(t(mat)/cs)*1e6 + 1)
  } else if (scale_type %in% c("tpm","fpkm","norm")) {
    mat <- log2(mat + 1)
  } else if (scale_type == "log2") {
    # already log
  } else stop("scale_type")
  keep <- rep(TRUE, nrow(mat))   # PATCH: no detection filter (paper convention)
  mat <- mat[keep, , drop=FALSE]
  sym <- rownames(mat)
  if (anyDuplicated(sym)) {
    o <- order(rowMeans(mat), decreasing=TRUE)
    mat <- mat[o, , drop=FALSE]; mat <- mat[!duplicated(rownames(mat)), , drop=FALSE]
  }
  mat
}

score_matrix <- function(expr, sets) {
  rk <- rankGenes(expr)
  sapply(sets, function(g) {
    gg <- intersect(g, rownames(expr))
    if (length(gg) < 5) return(rep(NA_real_, ncol(expr)))
    simpleScore(rk, upSet=gg, centerScore=TRUE)$TotalScore
  })
}

# --- all within-stratum label assignments (exhaustive if small, else sampled)
perm_assignments <- function(grp, strata, maxn=20000) {
  st <- unique(strata)
  lists <- lapply(st, function(s) {
    idx <- which(strata == s); g <- grp[idx]; nt <- sum(g == 1)
    cmb <- combn(length(idx), nt, simplify=FALSE)
    lapply(cmb, function(k) { v <- rep(0, length(idx)); v[k] <- 1; setNames(v, idx) })
  })
  total <- prod(sapply(lists, length))
  if (total <= maxn) {
    grid <- expand.grid(lapply(lists, seq_along))
    out <- apply(grid, 1, function(r) {
      v <- rep(NA, length(grp))
      for (j in seq_along(st)) { a <- lists[[j]][[r[j]]]; v[as.integer(names(a))] <- a }
      v })
    out <- t(out)
  } else {
    out <- t(replicate(maxn, { v <- rep(NA, length(grp))
      for (j in seq_along(st)) { a <- lists[[j]][[sample(length(lists[[j]]),1)]]; v[as.integer(names(a))] <- a }
      v }))
  }
  list(perms=out, space=total, exhaustive=(total <= maxn))
}

delta_of <- function(sc, grp, strata, paired_id=NULL) {
  if (!is.null(paired_id)) {
    d <- tapply(seq_along(sc), paired_id, function(i)
      mean(sc[i][grp[i]==1]) - mean(sc[i][grp[i]==0]))
    return(mean(d, na.rm=TRUE))
  }
  ds <- tapply(seq_along(sc), strata, function(i)
    mean(sc[i][grp[i]==1]) - mean(sc[i][grp[i]==0]))
  mean(ds, na.rm=TRUE)
}

score_dataset <- function(expr, grp, species, strata=NULL, paired_id=NULL,
                          nrand=1000, seed=1, maxn=20000) {
  set.seed(seed)
  if (is.null(strata)) strata <- rep(1, length(grp))
  sets <- get_sets(species)
  sc <- score_matrix(expr, sets)                     # samples x modules
  bg <- rownames(expr)
  pa <- perm_assignments(grp, strata, maxn)
  out <- list()
  for (m in MODULES) {
    s <- sc[, m]
    if (all(is.na(s))) { out[[m]] <- data.frame(module=m, n_genes_hit=NA); next }
    d <- delta_of(s, grp, strata, paired_id)
    # SE / CI
    if (!is.null(paired_id)) {
      dv <- tapply(seq_along(s), paired_id, function(i) mean(s[i][grp[i]==1])-mean(s[i][grp[i]==0]))
      se <- sd(dv)/sqrt(length(dv)); df <- length(dv)-1
    } else {
      n1 <- sum(grp==1); n0 <- sum(grp==0)
      v1 <- if (n1>1) var(s[grp==1]) else NA; v0 <- if (n0>1) var(s[grp==0]) else NA
      se <- suppressWarnings(sqrt(v1/n1 + v0/n0)); df <- n1+n0-2
    }
    ci <- if (is.finite(se) && df>0) d + c(-1,1)*qt(0.975, df)*se else c(NA,NA)
    # within-study permutation
    pd <- apply(pa$perms, 1, function(g) delta_of(s, g, strata, paired_id))
    pperm <- mean(abs(pd) >= abs(d) - 1e-12)
    # size-matched random gene sets
    k <- length(intersect(sets[[m]], bg))
    rnd <- replicate(nrand, {
      gg <- sample(bg, k)
      rs <- simpleScore(rankGenes(expr), upSet=gg, centerScore=TRUE)$TotalScore
      delta_of(rs, grp, strata, paired_id) })
    z <- (d - mean(rnd))/sd(rnd)
    prand <- mean(abs(rnd - mean(rnd)) >= abs(d - mean(rnd)))
    out[[m]] <- data.frame(module=m, n_genes_hit=k, delta=d, se=se,
      ci_lo=ci[1], ci_hi=ci[2], perm_space=pa$space, perm_exhaustive=pa$exhaustive,
      p_perm=pperm, z_rand=z, p_rand=prand)
  }
  res <- do.call(rbind, out)
  # DE / amplitude
  if (sum(grp==1)>1 && sum(grp==0)>1) {
    des <- if (length(unique(strata))>1) model.matrix(~factor(strata)+grp) else model.matrix(~grp)
    fit <- eBayes(lmFit(expr, des), trend=TRUE)
    tt <- topTable(fit, coef="grp", number=Inf, sort.by="none")
    res$n_DE <- sum(tt$adj.P.Val < 0.05); res$rms_lfc <- sqrt(mean(tt$logFC^2))
  } else {
    lfc <- rowMeans(expr[, grp==1, drop=FALSE]) - rowMeans(expr[, grp==0, drop=FALSE])
    res$n_DE <- NA_integer_; res$rms_lfc <- sqrt(mean(lfc^2))
  }
  res$delta_over_rms <- res$delta/res$rms_lfc
  attr(res, "scores") <- sc
  res
}

