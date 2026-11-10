suppressMessages({library(singscore);library(org.Mm.eg.db);library(AnnotationDbi)})
set.seed(1)
C<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
D<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/项目二_交付学生审阅/05_provenance"
GS<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
EF<-read.csv(file.path(D,"B1_effector_genesets.csv"),stringsAsFactors=FALSE)
cat("effector set columns:",paste(names(EF),collapse=","),"\n"); print(head(EF,3))
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
sp<-SP[SP$study=="S1a_GSE274336",]
M<-read.csv(file.path(C,"GSE274336_counts.csv.gz"),row.names=1,check.names=FALSE)
m<-as.matrix(M); mode(m)<-"numeric"; m[!is.finite(m)]<-0
cc<-intersect(sp$sample,colnames(m)); m<-m[,cc]
g<-ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0); st<-sp$stratum[match(cc,sp$sample)]
cs<-colSums(m); L<-log2(t(t(m)/cs)*1e6+1)
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(L)),keytype="ENSEMBL",column="SYMBOL",multiVals="first")))
ok<-!is.na(s)&s!=""; X<-L[ok,,drop=FALSE]; rownames(X)<-s[ok]
o<-order(rowMeans(X),decreasing=TRUE); X<-X[o,]; X<-X[!duplicated(rownames(X)),]
R<-rankGenes(X)
dd<-function(hit){ sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
  mean(sapply(unique(st),function(u){k<-st==u;mean(sc[k&g==1])-mean(sc[k&g==0])})) }
perm<-function(hit,B=20000){ obs<-dd(hit); cnt<-0
  for(i in 1:B){ gp<-g; for(u in unique(st)){k<-which(st==u); gp[k]<-sample(g[k])}
    sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
    v<-mean(sapply(unique(st),function(u){q<-st==u;mean(sc[q&gp==1])-mean(sc[q&gp==0])}))
    if(abs(v)>=abs(obs)) cnt<-cnt+1 }
  (cnt+1)/(B+1) }
cat(sprintf("\nbackground (no detection filter) = %d genes\n\n",nrow(X)))
cat(sprintf("%-26s %4s %10s %9s %9s\n","module","nhit","delta","z_rand","p_rand"))
sets<-list()
for(k in c("ISG_ISGF3_m","YAP_m")) sets[[k]]<-unique(strsplit(GS$genes[GS$set==k],";")[[1]])
cn<-names(EF)
for(i in seq_len(nrow(EF))){
  nmv<-as.character(EF[i,1]); gv<-as.character(EF[i,ncol(EF)])
  if(grepl(";",gv)) sets[[nmv]]<-unique(strsplit(gv,";")[[1]])
}
for(nm in names(sets)){ hit<-intersect(sets[[nm]],rownames(X)); if(length(hit)<5) next
  obs<-dd(hit); nul<-replicate(1000, dd(sample(rownames(X),length(hit))))
  cat(sprintf("%-26s %4d %+10.5f %+9.2f %9.3f\n",nm,length(hit),obs,(obs-mean(nul))/sd(nul),mean(abs(nul)>=abs(obs)))) }
