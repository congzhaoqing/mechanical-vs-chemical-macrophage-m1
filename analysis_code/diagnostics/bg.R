suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
CACHE<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)
isgf3<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
prep<-function(mat,st){mat<-as.matrix(mat);mode(mat)<-"numeric";mat[!is.finite(mat)]<-0
  if(st=="counts"){cs<-colSums(mat);cs[cs==0]<-1;mat<-log2(t(t(mat)/cs)*1e6+1)}
  keep<-rowSums(mat>0)>=max(2,floor(0.25*ncol(mat)));mat<-mat[keep,,drop=FALSE]
  if(anyDuplicated(rownames(mat))){o<-order(rowMeans(mat),decreasing=TRUE);mat<-mat[o,,drop=FALSE]
    mat<-mat[!duplicated(rownames(mat)),,drop=FALSE]};mat}
by_sym<-function(m,s){k<-!is.na(s)&s!="";m<-m[k,,drop=FALSE];rownames(m)<-s[k];m}
full<-read.csv(file.path(CACHE,"GSE274020_20240802_Raw_Counts.csv.gz"),row.names=1,check.names=FALSE)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(full)),
     keytype="ENTREZID",column="SYMBOL",multiVals="first")))
full<-by_sym(as.matrix(full),s)
cc<-grep("^WT-DMSO",colnames(full),value=TRUE); grp<-ifelse(grepl("-LPS-",cc),1,0)
d1<-{X<-prep(full[,cc],"counts");sc<-simpleScore(rankGenes(X),upSet=intersect(isgf3,rownames(X)),centerScore=TRUE)$TotalScore
     mean(sc[grp==1])-mean(sc[grp==0])}
d2<-{X<-prep(full,"counts");X<-X[,cc];sc<-simpleScore(rankGenes(X),upSet=intersect(isgf3,rownames(X)),centerScore=TRUE)$TotalScore
     mean(sc[grp==1])-mean(sc[grp==0])}
Xf<-prep(full,"counts");scf<-simpleScore(rankGenes(Xf),upSet=intersect(isgf3,rownames(Xf)),centerScore=TRUE)$TotalScore
names(scf)<-colnames(Xf); d3<-mean(scf[cc][grp==1])-mean(scf[cc][grp==0])
cat(sprintf("subset-then-prep  : %+.5f\nprep-full-then-subset(rank within subset): %+.5f\nrank on full matrix, then contrast: %+.5f\npublished: +0.12450\n",d1,d2,d3))
