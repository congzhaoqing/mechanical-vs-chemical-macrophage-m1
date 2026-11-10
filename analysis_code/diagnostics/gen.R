suppressMessages({library(singscore);library(org.Mm.eg.db);library(org.Hs.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
PB<-read.csv(file.path(C,"_pubs.csv"),stringsAsFactors=FALSE)
idtype<-function(r){ r<-sub("\\..*$","",r)
  if(mean(grepl("^ENS",r))>0.5) "ENSEMBL" else if(mean(grepl("^(NM_|XM_|NR_|XR_)",r))>0.5) "REFSEQ"
  else if(mean(grepl("^[0-9]+$",r))>0.8) "ENTREZID" else "SYMBOL" }
rd<-function(f){ sep<-if(grepl("tsv|txt",f)) "\t" else ","
  read.csv(file.path(C,f),row.names=1,check.names=FALSE,sep=sep) }
cat(sprintf("%-18s %-7s %6s %6s %9s %9s %9s\n","study","idtype","bg_filt","bg_all","published","filtered","unfiltered"))
for(i in seq_len(nrow(PB))){
  st<-PB$study[i]; sp<-SP[SP$study==st,]
  db<-if(PB$species[i]=="mouse") org.Mm.eg.db else org.Hs.eg.db
  key<-if(PB$species[i]=="mouse") "ISG_ISGF3_m" else "ISG_ISGF3_h"
  isg<-unique(strsplit(GS$genes[GS$set==key],";")[[1]])
  M<-rd(PB$file[i]); it<-idtype(rownames(M))
  if(it=="SYMBOL"){ s<-rownames(M) } else {
    s<-unname(suppressMessages(mapIds(db,keys=sub("\\..*$","",rownames(M)),keytype=it,column="SYMBOL",multiVals="first"))) }
  cc<-intersect(sp$sample,colnames(M))
  if(length(cc)<nrow(sp)) { cat(sprintf("%-18s SAMPLE MISMATCH %d/%d\n",st,length(cc),nrow(sp))); next }
  g<-ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0)
  m0<-as.matrix(M)[,cc,drop=FALSE]; mode(m0)<-"numeric"; m0[!is.finite(m0)]<-0
  L<-if(PB$scale[i]=="counts"){cs<-colSums(m0);cs[cs==0]<-1;log2(t(t(m0)/cs)*1e6+1)} else log2(m0+1)
  ok<-!is.na(s)&s!=""
  mk<-function(X,sy){ rownames(X)<-sy; o<-order(rowMeans(X),decreasing=TRUE); X<-X[o,,drop=FALSE]
                      X[!duplicated(rownames(X)),,drop=FALSE] }
  Xall<-mk(L[ok,,drop=FALSE],s[ok])
  keep<-rowSums(L>0)>=max(2,floor(0.25*ncol(L))); Xf<-mk(L[keep&ok,,drop=FALSE],s[keep&ok])
  d<-function(X){ h<-intersect(isg,rownames(X)); if(length(h)<5) return(NA)
    sc<-simpleScore(rankGenes(X),upSet=h,centerScore=TRUE)$TotalScore
    st2<-sp$stratum[match(cc,sp$sample)]
    if(length(unique(st2))>1) mean(sapply(unique(st2),function(u){k<-st2==u;mean(sc[k&g==1])-mean(sc[k&g==0])}))
    else mean(sc[g==1])-mean(sc[g==0]) }
  cat(sprintf("%-18s %-7s %6d %6d %+9.5f %+9.5f %+9.5f\n",st,it,nrow(Xf),nrow(Xall),PB$published[i],d(Xf),d(Xall)))
}
