suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
C<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
M<-read.csv(file.path(C,"GSE239377_collagen_RNAseq_TPM.txt.gz"),row.names=1,check.names=FALSE,sep="\t")
M<-M[,setdiff(colnames(M),"external_gene_name")]
m<-as.matrix(M); mode(m)<-"numeric"; m[!is.finite(m)]<-0; L<-log2(m+1)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(L)),keytype="ENSEMBL",column="SYMBOL",multiVals="first")))
ok<-!is.na(s)&s!=""; L<-L[ok,]; rownames(L)<-s[ok]
o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,]; L<-L[!duplicated(rownames(L)),]
SC<-simpleScore(rankGenes(L),upSet=intersect(isg,rownames(L)),centerScore=TRUE)$TotalScore
cn<-colnames(L); sub<-sub("_.*","",cn); il4<-grepl("IL-4",cn)
PUB<--0.008268; best<-list(d=99,tag="")
gr<-c("TC","7mg","2mg")
for(a in gr) for(b in gr){ if(a==b) next
 for(mode in c("untx","il4","both_pooled","both_strat")){
  k<-(sub %in% c(a,b)); if(mode=="untx") k<-k&!il4; if(mode=="il4") k<-k&il4
  if(sum(k)<4) next
  g<-ifelse(sub[k]==b,1,0); st<-if(mode=="both_strat") ifelse(il4[k],"i","u") else rep("a",sum(k))
  d<-mean(sapply(unique(st),function(u){q<-st==u; mean(SC[k][q&g==1])-mean(SC[k][q&g==0])}))
  tag<-sprintf("%s vs %s (%s)",a,b,mode)
  if(abs(d-PUB)<abs(best$d-PUB)) best<-list(d=d,tag=tag)
  if(abs(d-PUB)<0.0006) cat(sprintf("   %-28s n=%2d delta %+0.5f  <<< MATCH\n",tag,sum(k),d))
 }}
# collagen pooled vs TC, stratified by IL4
for(fl in c(FALSE,TRUE)){ g<-ifelse(sub!="TC",1,0); if(fl) g<-1-g
  st<-ifelse(il4,"i","u")
  d<-mean(sapply(unique(st),function(u){q<-st==u; mean(SC[q&g==1])-mean(SC[q&g==0])}))
  tag<-sprintf("TC vs collagen pooled, strat IL4%s",ifelse(fl," flipped",""))
  if(abs(d-PUB)<abs(best$d-PUB)) best<-list(d=d,tag=tag)
  if(abs(d-PUB)<0.0006) cat(sprintf("   %-28s n=%2d delta %+0.5f  <<< MATCH\n",tag,length(g),d)) }
cat(sprintf("\npublished %+0.5f   closest variant: %s  delta %+0.5f  (gap %.5f)\n",PUB,best$tag,best$d,abs(best$d-PUB)))
