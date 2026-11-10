suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
cat("gene sets available:\n"); print(data.frame(set=GS$set, n=sapply(strsplit(GS$genes,";"),function(x)length(unique(x)))))
full<-read.csv(file.path(C,"GSE274020_20240802_Raw_Counts.csv.gz"),row.names=1,check.names=FALSE)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(full)),keytype="ENTREZID",column="SYMBOL",multiVals="first")))
k<-!is.na(s)&s!=""; m<-as.matrix(full)[k,,drop=FALSE]; rownames(m)<-s[k]
cc<-grep("^WT-DMSO",colnames(m),value=TRUE); g<-ifelse(grepl("-LPS-",cc),1,0)
m0<-m[,cc]; cs<-colSums(m0); L<-log2(t(t(m0)/cs)*1e6+1)
X<-L[rowSums(L>0)>=max(2,floor(0.25*ncol(L))),,drop=FALSE]
o<-order(rowMeans(X),decreasing=TRUE); X<-X[o,]; X<-X[!duplicated(rownames(X)),]
cat(sprintf("\ncontrast: %d samples (%d LPS vs %d vehicle), background %d genes\n",length(cc),sum(g),sum(g==0),nrow(X)))
R<-rankGenes(X)
pub<-c(ISG_ISGF3_m=0.124549, ISG_IFNa_m=0.081324, ISG_IFNg_m=0.099041, YAP_m=0.008767)
cat(sprintf("\n%-14s %4s %10s %10s %8s\n","module","nhit","rebuilt","published","ratio"))
for(sn in GS$set){
  gl<-unique(strsplit(GS$genes[GS$set==sn],";")[[1]]); hit<-intersect(gl,rownames(X))
  if(length(hit)<5) next
  sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
  d<-mean(sc[g==1])-mean(sc[g==0])
  p<-if(sn %in% names(pub)) pub[[sn]] else NA
  cat(sprintf("%-14s %4d %10.5f %10s %8s\n",sn,length(hit),d,
      ifelse(is.na(p),"-",sprintf("%.5f",p)), ifelse(is.na(p),"-",sprintf("%.2f",d/p))))
}
