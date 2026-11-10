suppressMessages({library(singscore);library(limma);library(org.Mm.eg.db);library(AnnotationDbi)})
CACHE<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-rbind(read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)[,c("set","genes")],
          read.csv(file.path(CACHE,"genesets_effector.csv"),stringsAsFactors=FALSE)[,c("set","genes")])
getg<-function(sn) unique(strsplit(GS$genes[GS$set==sn],";")[[1]])
prep<-function(m,st){m<-as.matrix(m);mode(m)<-"numeric";m[!is.finite(m)]<-0
  if(st=="counts"){cs<-colSums(m);cs[cs==0]<-1;m<-log2(t(t(m)/cs)*1e6+1)}
  m<-m[rowSums(m>0)>=max(2,floor(0.25*ncol(m))),,drop=FALSE]
  if(anyDuplicated(rownames(m))){o<-order(rowMeans(m),decreasing=TRUE);m<-m[o,,drop=FALSE]
    m<-m[!duplicated(rownames(m)),,drop=FALSE]};m}
bysym<-function(m,s){k<-!is.na(s)&s!="";m<-m[k,,drop=FALSE];rownames(m)<-s[k];m}
m<-read.csv("raw2/S1a_GSE274336__raw.csv",row.names=1,check.names=FALSE)
md<-read.csv("raw2/S1a_GSE274336__meta.csv",colClasses="character")
s<-unname(suppressMessages(mapIds(org.Mm.eg.db,keys=sub("\\..*$","",rownames(m)),keytype="ENSEMBL",column="SYMBOL",multiVals="first")))
X<-prep(bysym(as.matrix(m),s),"counts")
g<-as.integer(md$grp); st<-md$stratum
cat(sprintf("GSE274336: %d genes x %d samples; strata %s\n",nrow(X),ncol(X),paste(table(st),collapse="/")))
dof<-function(sc) mean(tapply(seq_along(sc),st,function(i)mean(sc[i][g[i]==1])-mean(sc[i][g[i]==0])))
set.seed(42); NP<-20000
P<-t(replicate(NP,{v<-g; for(k in unique(st)){idx<-which(st==k); v[idx]<-sample(g[idx])}; v}))
rk<-rankGenes(X); bg<-rownames(X)
for(mod in c("ISG_ISGF3","EFF_TNFA_NFKB","EFF_INFLAMM","EFF_FIELD_PANEL","YAP")){
  hit<-intersect(getg(paste0(mod,"_m")),bg); if(length(hit)<5) next
  sc<-simpleScore(rk,upSet=hit,centerScore=TRUE)$TotalScore; d<-dof(sc)
  pd<-apply(P,1,function(gg){v<-gg; mean(tapply(seq_along(sc),st,function(i)mean(sc[i][v[i]==1])-mean(sc[i][v[i]==0])))})
  pp<-mean(abs(pd)>=abs(d)-1e-12)
  rnd<-replicate(1000,{gg<-sample(bg,length(hit))
    dof(simpleScore(rk,upSet=gg,centerScore=TRUE)$TotalScore)})
  z<-(d-mean(rnd))/sd(rnd); pr<-mean(abs(rnd-mean(rnd))>=abs(d-mean(rnd)))
  cat(sprintf("  %-16s n=%3d  delta %+0.5f  p_perm %.4f  z_rand %+5.2f  p_rand %.3f\n",
      mod,length(hit),d,pp,z,pr))
}
des<-model.matrix(~factor(st)+g); fit<-eBayes(lmFit(X,des),trend=TRUE)
tt<-topTable(fit,coef="g",number=Inf,sort.by="none")
cat(sprintf("\nDE genes FDR<0.05: %d  (published 3,101)\n",sum(tt$adj.P.Val<0.05)))
