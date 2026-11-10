suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
set.seed(1)
C<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
sp<-SP[SP$study=="S1a_GSE274336",]
M<-read.csv(file.path(C,"GSE274336_counts.csv.gz"),row.names=1,check.names=FALSE)
m<-as.matrix(M); mode(m)<-"numeric"; m[!is.finite(m)]<-0
cc<-intersect(sp$sample,colnames(m)); m<-m[,cc]
g<-ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0); st<-sp$stratum[match(cc,sp$sample)]
cs<-colSums(m); L<-log2(t(t(m)/cs)*1e6+1)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(L)),keytype="ENSEMBL",column="SYMBOL",multiVals="first")))
ok<-!is.na(s)&s!=""; X<-L[ok,,drop=FALSE]; rownames(X)<-s[ok]
o<-order(rowMeans(X),decreasing=TRUE); X<-X[o,]; X<-X[!duplicated(rownames(X)),]
R<-rankGenes(X)
dd<-function(hit){ sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
  mean(sapply(unique(st),function(u){k<-st==u;mean(sc[k&g==1])-mean(sc[k&g==0])})) }
for(key in c("ISG_ISGF3_m","YAP_m")){
  isg<-unique(strsplit(GS$genes[GS$set==key],";")[[1]]); hit<-intersect(isg,rownames(X))
  obs<-dd(hit); nul<-replicate(1000, dd(sample(rownames(X),length(hit))))
  cat(sprintf("%-12s bg=%d nhit=%3d  delta %+0.5f  z_rand %+0.2f  p_rand %.3f\n",
    key,nrow(X),length(hit),obs,(obs-mean(nul))/sd(nul),mean(abs(nul)>=abs(obs))))
}
cat("\noriginal reported for GSE274336: delta -0.0025, z -0.98, p_random 0.52\n")
