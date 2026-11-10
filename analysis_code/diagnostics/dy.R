suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
full<-read.csv(file.path(C,"GSE274020_20240802_Raw_Counts.csv.gz"),row.names=1,check.names=FALSE)
cat(sprintf("raw matrix rows (Entrez-indexed): %d\n",nrow(full)))
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(full)),keytype="ENTREZID",column="SYMBOL",multiVals="first")))
cat(sprintf("rows with a symbol: %d   unique symbols: %d\n",sum(!is.na(s)&s!=""),length(unique(s[!is.na(s)&s!=""]))))
cc<-grep("^WT-DMSO",colnames(full),value=TRUE); g<-ifelse(grepl("-LPS-",cc),1,0)
sub0<-as.matrix(full)[,cc]; cs<-colSums(sub0); Lall<-log2(t(t(sub0)/cs)*1e6+1)
D<-function(X,lab){
  hit<-intersect(isg,rownames(X))
  sc<-simpleScore(rankGenes(X),upSet=hit,centerScore=TRUE)$TotalScore
  d<-mean(sc[g==1])-mean(sc[g==0])
  cat(sprintf("  %-46s bg=%6d  nhit=%3d  delta %+0.5f  ratio_to_pub %.2f\n",lab,nrow(X),length(hit),d,d/0.124549))
}
cat("\npublished delta = +0.12455   (predicted background if 1/N scaling: ~25,900)\n\n")
# 1. current convention
X<-Lall[rowSums(Lall>0)>=2,,drop=FALSE]; rn<-s[match(rownames(X),rownames(full))]
X1<-X[!is.na(rn)&rn!="",,drop=FALSE]; r1<-rn[!is.na(rn)&rn!=""]; rownames(X1)<-r1
o<-order(rowMeans(X1),decreasing=TRUE); X1<-X1[o,]; X1<-X1[!duplicated(rownames(X1)),]; D(X1,"symbol-mapped, detected>=2, dedup  [current]")
# 2. keep ALL rows incl all-zero, symbol-mapped, dedup
Xa<-Lall; rn<-s[match(rownames(Xa),rownames(full))]
X2<-Xa[!is.na(rn)&rn!="",,drop=FALSE]; r2<-rn[!is.na(rn)&rn!=""]; rownames(X2)<-r2
o<-order(rowMeans(X2),decreasing=TRUE); X2<-X2[o,]; X2<-X2[!duplicated(rownames(X2)),]; D(X2,"symbol-mapped, NO detection filter, dedup")
# 3. keep all rows, unmapped retained under original id, no dedup
rn3<-s[match(rownames(Lall),rownames(full))]; rn3[is.na(rn3)|rn3==""]<-rownames(Lall)[is.na(rn3)|rn3==""]
X3<-Lall; rownames(X3)<-make.unique(rn3); D(X3,"ALL rows, unmapped kept, no dedup")
# 4. all rows, unmapped kept, detected>=2
X4<-X3[rowSums(X3>0)>=2,,drop=FALSE]; D(X4,"ALL rows, unmapped kept, detected>=2")
