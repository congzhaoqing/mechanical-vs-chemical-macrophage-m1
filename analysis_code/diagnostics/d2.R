suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
CACHE<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS<-read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
bysym<-function(m,s){k<-!is.na(s)&s!="";m<-m[k,,drop=FALSE];rownames(m)<-s[k];m}
prep<-function(m,st){m<-as.matrix(m);mode(m)<-"numeric";m[!is.finite(m)]<-0
  if(st=="counts"){cs<-colSums(m);cs[cs==0]<-1;m<-log2(t(t(m)/cs)*1e6+1)
  } else if(st %in% c("tpm","fpkm","norm")) m<-log2(m+1)
  m<-m[rowSums(m>0)>=max(2,floor(0.25*ncol(m))),,drop=FALSE]
  if(anyDuplicated(rownames(m))){o<-order(rowMeans(m),decreasing=TRUE);m<-m[o,,drop=FALSE]
    m<-m[!duplicated(rownames(m)),,drop=FALSE]};m}
D<-function(X,grp){sc<-simpleScore(rankGenes(X),upSet=intersect(isg,rownames(X)),centerScore=TRUE)$TotalScore
  mean(sc[grp==1])-mean(sc[grp==0])}
run<-function(file,key,kt,pat,lpspat,pub){
  full<-read.csv(file.path(CACHE,file),row.names=1,check.names=FALSE,sep=if(grepl("txt",file))"\t" else ",")
  ids<-sub("\\..*$","",rownames(full))
  s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=ids,keytype=kt,column="SYMBOL",multiVals="first")))
  full<-bysym(as.matrix(full),s); cc<-grep(pat,colnames(full),value=TRUE)
  g<-ifelse(grepl(lpspat,cc),1,0)
  cat(sprintf("\n%s  (published %+.4f)\n",key,pub))
  for(st in c("counts","norm")) cat(sprintf("   scale=%-7s delta %+0.5f\n",st,D(prep(full[,cc],st),g)))
}
run("GSE274020_20240802_Raw_Counts.csv.gz","GSE274020 WT-DMSO LPS 6h","ENTREZID","^WT-DMSO","-LPS-",0.1245)
run("GSE253423_raw_read_quant.table.txt.gz","GSE253423 WT LPS 4h","REFSEQ","^WT._(0h|4h)$","_4h$",0.1380)
