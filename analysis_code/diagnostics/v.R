suppressMessages({library(singscore);library(org.Mm.eg.db);library(org.Hs.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
prep<-function(M,sc,db,it){
  m<-as.matrix(M); mode(m)<-"numeric"; m[!is.finite(m)]<-0
  L<-if(sc=="counts"){cs<-colSums(m);cs[cs==0]<-1;log2(t(t(m)/cs)*1e6+1)} else log2(m+1)
  s<-if(it=="SYMBOL") rownames(L) else unname(suppressMessages(mapIds(db,keys=sub("\\..*$","",rownames(L)),keytype=it,column="SYMBOL",multiVals="first")))
  ok<-!is.na(s)&s!=""; L<-L[ok,,drop=FALSE]; rownames(L)<-s[ok]
  o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,,drop=FALSE]; L[!duplicated(rownames(L)),,drop=FALSE] }
sco<-function(X,isg) simpleScore(rankGenes(X),upSet=intersect(isg,rownames(X)),centerScore=TRUE)$TotalScore
dl<-function(sc,g,st){ if(length(unique(st))>1) mean(sapply(unique(st),function(u){k<-st==u;mean(sc[k&g==1])-mean(sc[k&g==0])}))
                       else mean(sc[g==1])-mean(sc[g==0]) }
show<-function(tag,pub,sc,g,st) cat(sprintf("   %-34s n=%2d  delta %+0.5f   %s\n",tag,length(g),dl(sc,g,st),
       ifelse(abs(dl(sc,g,st)-pub)<0.0005,"<<< MATCH","")))

cat("=== S4_GSE273834  published -0.08406 ===\n")
isg<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_m"],";")[[1]])
M<-read.csv(file.path(C,"GSE273834_Salmon.norm.tpm.csv.gz"),row.names=1,check.names=FALSE)
X<-prep(M,"tpm",org.Mm.eg.db,"SYMBOL"); s<-sco(X,isg); names(s)<-colnames(X)
cn<-colnames(X); st<-ifelse(grepl("^WT",cn),"WT","KO"); il4<-grepl("IL4",cn); hard<-grepl("_St_",cn)
k<-!il4; show("St vs So, strat genotype [current]",-0.08406,s[k],as.integer(hard[k]),st[k])
k<-!il4 & st=="WT"; show("WT only, St vs So",-0.08406,s[k],as.integer(hard[k]),rep("a",sum(k)))
k<-!il4 & st=="KO"; show("KO only, St vs So",-0.08406,s[k],as.integer(hard[k]),rep("a",sum(k)))
k<-!il4; show("So vs St (flipped)",-0.08406,s[k],as.integer(!hard[k]),st[k])
k<-rep(TRUE,length(cn)); show("all incl IL4, St vs So",-0.08406,s[k],as.integer(hard[k]),st[k])
k<-st=="WT"; show("WT, IL4 vs untreated",-0.08406,s[k],as.integer(il4[k]),rep("a",sum(k)))
k<-rep(TRUE,length(cn)); show("IL4 vs untreated, strat genotype",-0.08406,s[k],as.integer(il4[k]),st[k])

cat("\n=== S5_GSE239377  published -0.00827 ===\n")
M<-read.csv(file.path(C,"GSE239377_collagen_RNAseq_TPM.txt.gz"),row.names=1,check.names=FALSE,sep="\t")
M<-M[,setdiff(colnames(M),"external_gene_name")]
X<-prep(M,"tpm",org.Mm.eg.db,"ENSEMBL"); s<-sco(X,isg); cn<-colnames(X)
tc<-grepl("^TC",cn); m7<-grepl("^7mg",cn); m2<-grepl("^2mg",cn); il4<-grepl("IL-4",cn)
A<-function(sel,gv,tag) show(tag,-0.00827,s[sel],as.integer(gv[sel]),rep("a",sum(sel)))
A(!il4, m7|m2, "TC vs 7mg+2mg, untx [current]")
A(!il4, tc, "7mg+2mg vs TC (flipped)")
A(!il4 & (tc|m7), m7, "TC vs 7mg only")
A(!il4 & (tc|m2), m2, "TC vs 2mg only")
A(!il4 & (m7|m2), m7, "2mg vs 7mg")
A(rep(TRUE,length(cn)), m7|m2, "all incl IL4, TC vs collagen")
A(il4, m7|m2, "IL-4 arm only, TC vs collagen")
A(tc, il4, "TC only, IL-4 vs untx")
A(rep(TRUE,length(cn)), il4, "IL-4 vs untx, all substrates")

cat("\n=== S2_GSE321887  published -0.01143 ===\n")
isgh<-unique(strsplit(GS$genes[GS$set=="ISG_ISGF3_h"],";")[[1]])
M<-read.csv(file.path(C,"GSE321887_RNAseq_countMatrix.tsv.gz"),row.names=1,check.names=FALSE,sep="\t")
X<-prep(M,"counts",org.Hs.eg.db,"ENSEMBL"); s<-sco(X,isgh)
idx<-as.integer(sub(".*Sample_","",colnames(X))); names(s)<-idx
don<-((idx-1)%/%3)+1; pos<-((idx-1)%%3)+1
for(cmp in list(c(1,3),c(1,2),c(2,3))) for(fl in c(FALSE,TRUE)){
  k<-pos %in% cmp; gv<-ifelse(pos==cmp[2],1,0); if(fl) gv<-1-gv
  show(sprintf("pos %d vs %d%s",cmp[1],cmp[2],ifelse(fl," flipped","")),-0.01143,s[k],gv[k],paste0("D",don[k]))
}
