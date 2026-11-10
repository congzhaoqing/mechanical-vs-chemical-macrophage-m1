suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
go<-function(file,sep,it,study){
  M<-read.csv(file.path(C,file),row.names=1,check.names=FALSE,sep=sep)
  m<-as.matrix(M); mode(m)<-"numeric"; m[!is.finite(m)]<-0
  sp<-SP[SP$study==study,]; cc<-intersect(sp$sample,colnames(m))
  if(length(cc)<nrow(sp)){cat(sprintf("%s: matched %d/%d samples\n",study,length(cc),nrow(sp)));return(NULL)}
  m<-m[,cc,drop=FALSE]; cs<-colSums(m); cs[cs==0]<-1; L<-log2(t(t(m)/cs)*1e6+1)
  s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(L)),keytype=it,column="SYMBOL",multiVals="first")))
  ok<-!is.na(s)&s!=""; L<-L[ok,]; rownames(L)<-s[ok]
  o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,]; L<-L[!duplicated(rownames(L)),]
  sc<-simpleScore(rankGenes(L),upSet=intersect(isg,rownames(L)),centerScore=TRUE)$TotalScore
  g<-ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0); st<-sp$stratum[match(cc,sp$sample)]
  d<-if(length(unique(st))>1) mean(sapply(unique(st),function(u){k<-st==u;mean(sc[k&g==1])-mean(sc[k&g==0])})) else mean(sc[g==1])-mean(sc[g==0])
  v<-(var(sc[g==1])/sum(g==1))+(var(sc[g==0])/sum(g==0))
  cat(sprintf("%-16s bg=%6d n=%2d  delta %+0.5f  se %0.5f\n",study,nrow(L),length(cc),d,sqrt(v)))
  c(d=d,se=sqrt(v)) }
a<-go("GSE274336_counts.csv.gz",",","SYMBOL","S1a_GSE274336")
b<-go("GSE274019_20240603_counts.txt.gz","\t","SYMBOL","S1b_GSE274019")
if(!is.null(a)&&!is.null(b)){
 cat(sprintf("\nS1 merge options (published +0.004305):\n  simple mean   %+0.5f\n  inv-variance  %+0.5f\n",
   mean(c(a[["d"]],b[["d"]])), (a[["d"]]/a[["se"]]^2+b[["d"]]/b[["se"]]^2)/(1/a[["se"]]^2+1/b[["se"]]^2))) }
