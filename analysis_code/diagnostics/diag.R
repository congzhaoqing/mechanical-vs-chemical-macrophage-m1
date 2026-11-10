suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
CACHE<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS<-read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
prep<-function(m,st){m<-as.matrix(m);mode(m)<-"numeric";m[!is.finite(m)]<-0
  if(st=="counts"){cs<-colSums(m);cs[cs==0]<-1;m<-log2(t(t(m)/cs)*1e6+1)}
  m<-m[rowSums(m>0)>=max(2,floor(0.25*ncol(m))),,drop=FALSE]
  if(anyDuplicated(rownames(m))){o<-order(rowMeans(m),decreasing=TRUE);m<-m[o,,drop=FALSE]
    m<-m[!duplicated(rownames(m)),,drop=FALSE]};m}
bysym<-function(m,s){k<-!is.na(s)&s!="";m<-m[k,,drop=FALSE];rownames(m)<-s[k];m}
D<-function(X,grp){sc<-simpleScore(rankGenes(X),upSet=intersect(isg,rownames(X)),centerScore=TRUE)$TotalScore
  mean(sc[grp==1])-mean(sc[grp==0])}
full<-read.csv(file.path(CACHE,"GSE274020_20240802_Raw_Counts.csv.gz"),row.names=1,check.names=FALSE)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(full)),keytype="ENTREZID",column="SYMBOL",multiVals="first")))
full<-bysym(as.matrix(full),s); cn<-colnames(full)
cat("GSE274020 — published anchor +0.12450\n")
variants<-list(
 "WT-DMSO LPS vs veh"      = cn[grepl("^WT-DMSO",cn)],
 "WT all arms LPS vs none" = cn[grepl("^WT-",cn)],
 "all genotypes LPS vs none"= cn,
 "KO-DMSO LPS vs veh"      = cn[grepl("^NRF2KO-DMSO",cn)],
 "WT+KO DMSO only"         = cn[grepl("DMSO",cn)])
for(nm in names(variants)){ cc<-variants[[nm]]; g<-ifelse(grepl("-LPS-",cc),1,0)
  X<-prep(full[,cc],"counts"); cat(sprintf("  %-28s n=%2d  delta %+0.5f\n",nm,length(cc),D(X,g))) }
