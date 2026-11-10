suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
CACHE<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
bysym<-function(m,s){k<-!is.na(s)&s!="";m<-m[k,,drop=FALSE];rownames(m)<-s[k];m}
full<-read.csv(file.path(CACHE,"GSE274020_20240802_Raw_Counts.csv.gz"),row.names=1,check.names=FALSE)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(full)),keytype="ENTREZID",column="SYMBOL",multiVals="first")))
fm<-bysym(as.matrix(full),s); cc<-grep("^WT-DMSO",colnames(fm),value=TRUE); g<-ifelse(grepl("-LPS-",cc),1,0)
m0<-fm[,cc]; cs<-colSums(m0); L<-log2(t(t(m0)/cs)*1e6+1)
cat("background-size sensitivity of the SAME contrast (GSE274020 WT-DMSO LPS 6h)\n")
cat("published value +0.1245\n\n")
for(frac in c(0,0.10,0.25,0.50,0.75)){
  thr<-max(2,floor(frac*ncol(L))); X<-L[rowSums(L>0)>=thr,,drop=FALSE]
  if(anyDuplicated(rownames(X))){o<-order(rowMeans(X),decreasing=TRUE);X<-X[o,];X<-X[!duplicated(rownames(X)),]}
  sc<-simpleScore(rankGenes(X),upSet=intersect(isg,rownames(X)),centerScore=TRUE)$TotalScore
  cat(sprintf("  detect-in >= %4.0f%% of samples -> background %6d genes, delta %+0.5f\n",
      frac*100, nrow(X), mean(sc[g==1])-mean(sc[g==0])))
}
# unfiltered, INCLUDING genes that failed symbol mapping (kept as their original id)
m2<-as.matrix(full)[,cc]; L2<-log2(t(t(m2)/colSums(m2))*1e6+1)
X2<-L2[rowSums(L2>0)>=2,,drop=FALSE]
rn<-s[match(rownames(X2),rownames(full))]; rn[is.na(rn)]<-rownames(X2)[is.na(rn)]
rownames(X2)<-rn; X2<-X2[!duplicated(rn),]
sc<-simpleScore(rankGenes(X2),upSet=intersect(isg,rownames(X2)),centerScore=TRUE)$TotalScore
cat(sprintf("\n  keeping unmapped ids in background -> %6d genes, delta %+0.5f\n",nrow(X2),mean(sc[g==1])-mean(sc[g==0])))
