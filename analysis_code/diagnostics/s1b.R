suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
C<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
score<-function(L,g,st){ sc<-simpleScore(rankGenes(L),upSet=intersect(isg,rownames(L)),centerScore=TRUE)$TotalScore
  d<-if(length(unique(st))>1) mean(sapply(unique(st),function(u){k<-st==u;mean(sc[k&g==1])-mean(sc[k&g==0])})) else mean(sc[g==1])-mean(sc[g==0])
  list(d=d, se=sqrt(var(sc[g==1])/sum(g==1)+var(sc[g==0])/sum(g==0)), n=length(g), bg=nrow(L)) }
mk<-function(m,it){ cs<-colSums(m); cs[cs==0]<-1; L<-log2(t(t(m)/cs)*1e6+1)
  s<-if(it=="SYMBOL") rownames(L) else unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(L)),keytype=it,column="SYMBOL",multiVals="first")))
  ok<-!is.na(s)&s!=""; L<-L[ok,,drop=FALSE]; rownames(L)<-s[ok]
  o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,,drop=FALSE]; L[!duplicated(rownames(L)),,drop=FALSE] }
go<-function(file,sep,it,study){
  M<-read.csv(file.path(C,file),row.names=1,check.names=FALSE,sep=sep)
  m<-as.matrix(M); mode(m)<-"numeric"; m[!is.finite(m)]<-0
  sp<-SP[SP$study==study,]; cc<-intersect(sp$sample,colnames(m))
  if(length(cc)<nrow(sp)){cat(sprintf("%-16s matched %d/%d\n",study,length(cc),nrow(sp)));return(NULL)}
  r<-score(mk(m[,cc,drop=FALSE],it), ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0), sp$stratum[match(cc,sp$sample)])
  cat(sprintf("%-16s bg=%6d n=%2d  delta %+0.5f  se %0.5f\n",study,r$bg,r$n,r$d,r$se)); r }
a<-go("GSE274336_counts.csv.gz",",","SYMBOL","S1a_GSE274336")
b<-go("GSE274019_20240603_counts.txt.gz","\t","SYMBOL","S1b_GSE274019")
if(!is.null(a)&&!is.null(b)) cat(sprintf("\nS1 published +0.004305 |  mean %+0.5f  |  inv-var %+0.5f\n",
  mean(c(a$d,b$d)), (a$d/a$se^2+b$d/b$se^2)/(1/a$se^2+1/b$se^2)))
# S6 GSE291556 from RAW.tar
td<-tempfile(); dir.create(td); untar(file.path(C,"GSE291556_RAW.tar"),exdir=td)
fs<-list.files(td,full.names=TRUE); nm<-sub("^GSM[0-9]+_","",sub("\\.txt\\.gz$","",basename(fs)))
ms<-lapply(fs,function(f){d<-read.csv(f,sep="\t",header=TRUE,row.names=1,check.names=FALSE); d[,ncol(d),drop=TRUE]})
gg<-read.csv(fs[1],sep="\t",header=TRUE,row.names=1,check.names=FALSE)
m6<-do.call(cbind,ms); rownames(m6)<-rownames(gg); colnames(m6)<-nm
cat(sprintf("\nGSE291556 matrix %d x %d : %s\n",nrow(m6),ncol(m6),paste(nm,collapse=", ")))
k<-nm %in% c("low_1","low_2","high_1","high_2")
r<-score(mk(m6[,k,drop=FALSE],"SYMBOL"), as.integer(grepl("^high",nm[k])), rep("a",sum(k)))
cat(sprintf("S6 high vs low   bg=%6d n=%2d  delta %+0.5f   (published -0.005883)\n",r$bg,r$n,r$d))
