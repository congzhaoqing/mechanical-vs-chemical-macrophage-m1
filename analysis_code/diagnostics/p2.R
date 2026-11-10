suppressMessages({library(singscore);library(org.Mm.eg.db);library(org.Hs.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
G1<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
G2<-read.csv(file.path(C,"genesets_effector.csv"),stringsAsFactors=FALSE); names(G2)[1]<-"set"
GS<-rbind(G1[,c("set","genes")],G2[,c("set","genes")])
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
MODS<-c(ISG="ISG_ISGF3",EFF_TNFA="EFF_TNFA_NFKB",EFF_INFL="EFF_INFLAMM",EFF_PANEL="EFF_FIELD_PANEL")
mk<-function(m,it,sc,sp){ mode(m)<-"numeric"; m[!is.finite(m)]<-0
  L<-if(sc=="counts"){cs<-colSums(m);cs[cs==0]<-1;log2(t(t(m)/cs)*1e6+1)} else if(sc=="log2") m else log2(m+1)
  db<-if(sp=="human") org.Hs.eg.db else org.Mm.eg.db
  s<-if(it=="SYMBOL") rownames(L) else unname(suppressMessages(mapIds(db,keys=sub("\\..*$","",rownames(L)),keytype=it,column="SYMBOL",multiVals="first")))
  ok<-!is.na(s)&s!=""; L<-L[ok,,drop=FALSE]; rownames(L)<-s[ok]
  o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,,drop=FALSE]; L[!duplicated(rownames(L)),,drop=FALSE] }
emit<-function(st,X,g,sv,sp){ R<-rankGenes(X); sfx<-if(sp=="human")"_h" else "_m"; v<-c()
  for(nm in names(MODS)){ gl<-unique(strsplit(GS$genes[GS$set==paste0(MODS[[nm]],sfx)],";")[[1]])
    hit<-intersect(gl,rownames(X)); if(length(hit)<5){v[nm]<-NA;next}
    sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
    v[nm]<-if(length(unique(sv))>1) mean(sapply(unique(sv),function(u){k<-sv==u;mean(sc[k&g==1])-mean(sc[k&g==0])})) else mean(sc[g==1])-mean(sc[g==0]) }
  cat(sprintf("%-18s bg=%6d n=%2d  ISG %+0.5f  TNFA %+0.5f  INFL %+0.5f  PANEL %+0.5f\n",st,nrow(X),length(g),v[1],v[2],v[3],v[4])); v }
res<-list()
# GSE112372 human LPS 3h anchor
M<-read.csv(file.path(C,"GSE112372_raw_gene_counts.txt.gz"),row.names=1,check.names=FALSE,sep="\t")
sp<-SP[SP$study=="A_GSE112372_LPS3h",]; cc<-intersect(sp$sample,colnames(M))
cat(sprintf("GSE112372 matched %d/%d\n",length(cc),nrow(sp)))
if(length(cc)==nrow(sp)) res$A_GSE112372_LPS3h<-emit("A_GSE112372_LPS3h",mk(as.matrix(M)[,cc,drop=FALSE],"ENSEMBL","counts","human"),
  ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0),sp$stratum[match(cc,sp$sample)],"human")
# GSE312997
M<-read.csv(file.path(C,"GSE312997_sheet1_counts.csv"),row.names=1,check.names=FALSE)
sp<-SP[SP$study=="S11_GSE312997",]; cc<-intersect(sp$sample,colnames(M))
cat(sprintf("GSE312997 matched %d/%d\n",length(cc),nrow(sp)))
if(length(cc)==nrow(sp)) res$S11_GSE312997<-emit("S11_GSE312997",mk(as.matrix(M)[,cc,drop=FALSE],"SYMBOL","counts","human"),
  ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0),sp$stratum[match(cc,sp$sample)],"human")
# GSE291556 from extracted tar, map GSM -> file
fs<-sort(list.files("/tmp/g291",full.names=TRUE)); gsm<-sub("_.*","",basename(fs))
ms<-lapply(fs,function(f) read.csv(f,sep="\t",header=TRUE,row.names=1,check.names=FALSE)[,1])
m6<-do.call(cbind,ms); rownames(m6)<-rownames(read.csv(fs[1],sep="\t",header=TRUE,row.names=1,check.names=FALSE)); colnames(m6)<-gsm
sp<-SP[SP$study=="S6_GSE291556",]; cc<-intersect(sp$sample,colnames(m6))
cat(sprintf("GSE291556 matched %d/%d\n",length(cc),nrow(sp)))
if(length(cc)==nrow(sp)) res$S6_GSE291556<-emit("S6_GSE291556",mk(m6[,cc,drop=FALSE],"SYMBOL","norm","mouse"),
  ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0),rep("a",length(cc)),"mouse")
saveRDS(res,"/tmp/p2.rds")
