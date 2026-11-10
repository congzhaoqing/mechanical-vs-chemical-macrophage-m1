suppressMessages({library(limma);library(org.Mm.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
sp<-SP[SP$study=="S1a_GSE274336",]
M<-read.csv(file.path(C,"GSE274336_counts.csv.gz"),row.names=1,check.names=FALSE)
m<-as.matrix(M); mode(m)<-"numeric"; m[!is.finite(m)]<-0
cc<-intersect(sp$sample,colnames(m)); m<-m[,cc]
g<-ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0); st<-sp$stratum[match(cc,sp$sample)]
cs<-colSums(m); L<-log2(t(t(m)/cs)*1e6+1)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(L)),keytype="ENSEMBL",column="SYMBOL",multiVals="first")))
ok<-!is.na(s)&s!=""
run<-function(X,tag){ o<-order(rowMeans(X),decreasing=TRUE); X<-X[o,,drop=FALSE]; X<-X[!duplicated(rownames(X)),,drop=FALSE]
  des<-model.matrix(~factor(st)+g); fit<-eBayes(lmFit(X,des),trend=TRUE)
  tt<-topTable(fit,coef="g",number=Inf,sort.by="none")
  cat(sprintf("  %-22s background %6d   DE at FDR<0.05 = %5d\n",tag,nrow(X),sum(tt$adj.P.Val<0.05))) }
X1<-L[ok,,drop=FALSE]; rownames(X1)<-s[ok]; run(X1,"NO detection filter")
keep<-rowSums(L>0)>=max(2,floor(0.25*ncol(L))); X2<-L[keep&ok,,drop=FALSE]; rownames(X2)<-s[keep&ok]; run(X2,"with detection filter")
cat("\n  original reported: 3,101   |   my B1 step-8 value: 4,318\n")
