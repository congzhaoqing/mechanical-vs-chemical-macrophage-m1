suppressMessages({library(singscore);library(org.Mm.eg.db);library(org.Hs.eg.db);library(AnnotationDbi)})
CACHE <- "/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS <- rbind(read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)[,c("set","genes")],
            read.csv(file.path(CACHE,"genesets_effector.csv"),stringsAsFactors=FALSE)[,c("set","genes")])
getg <- function(sn) unique(strsplit(GS$genes[GS$set==sn],";")[[1]])
MODS <- c("ISG_ISGF3","ISG_IFNA","ISG_IFNG","YAP","EFF_TNFA_NFKB","EFF_INFLAMM","EFF_FIELD_PANEL")
prep_expr <- function(mat, st){ mat<-as.matrix(mat); mode(mat)<-"numeric"; mat[!is.finite(mat)]<-0
  if(st=="counts"){cs<-colSums(mat);cs[cs==0]<-1;mat<-log2(t(t(mat)/cs)*1e6+1)
  } else if(st %in% c("tpm","fpkm","norm")) mat<-log2(mat+1)
  keep<-rowSums(mat>0)>=max(2,floor(0.25*ncol(mat))); mat<-mat[keep,,drop=FALSE]
  if(anyDuplicated(rownames(mat))){o<-order(rowMeans(mat),decreasing=TRUE);mat<-mat[o,,drop=FALSE]
    mat<-mat[!duplicated(rownames(mat)),,drop=FALSE]}; mat }
by_sym<-function(m,s){k<-!is.na(s)&s!="";m<-m[k,,drop=FALSE];rownames(m)<-s[k];m}
KT <- list(S1a_GSE274336=c("ENSEMBL","mouse"),S2_GSE321887=c("ENSEMBL","human"),
  S7_GSE206771=c("ENSEMBL","mouse"),S8_GSE301422=c("ENSEMBL","human"),
  S1b_GSE274019=c("ENTREZID","mouse"),A_GSE274020_LPS6h=c("ENTREZID","mouse"),
  A_GSE253423_LPS4h=c("REFSEQ","mouse"),A_GSE253423_LPS18h=c("REFSEQ","mouse"),
  A_GSE99445_IFNg=c("ENTREZID","human"),A_GSE99446_IFNg=c("ENTREZID","human"),
  A_GSE99445_porousTi=c("ENTREZID","human"))
SCALE <- list(A_GSE99445_IFNg="log2",A_GSE99446_IFNg="log2",A_GSE99445_porousTi="log2")
res<-list()
for(dirn in c("raw2","raw3")) for(f in list.files(dirn,pattern="__raw.csv$")){
  st<-sub("__raw.csv","",f); kt<-KT[[st]]; if(is.null(kt)) next
  m<-read.csv(file.path(dirn,f),row.names=1,check.names=FALSE)
  md<-read.csv(file.path(dirn,paste0(st,"__meta.csv")),colClasses="character")
  ids<-sub("\\..*$","",rownames(m)); if(kt[1]=="ENTREZID") ids<-sub("_at$","",ids)
  db<-if(kt[2]=="mouse") org.Mm.eg.db else org.Hs.eg.db
  s<-unname(suppressMessages(mapIds(db,keys=ids,keytype=kt[1],column="SYMBOL",multiVals="first")))
  X<-prep_expr(by_sym(as.matrix(m),s), if(!is.null(SCALE[[st]])) SCALE[[st]] else "counts")
  grp<-as.integer(md$grp); strata<-md$stratum; paired<-md$paired
  rk<-rankGenes(X); sfx<-if(kt[2]=="human")"_h" else "_m"
  for(mm in MODS){ hit<-intersect(getg(paste0(mm,sfx)),rownames(X)); if(length(hit)<5) next
    sc<-simpleScore(rk,upSet=hit,centerScore=TRUE)$TotalScore
    if(any(paired!="")){dv<-tapply(seq_along(sc),paired,function(i)mean(sc[i][grp[i]==1])-mean(sc[i][grp[i]==0]))
      d<-mean(dv);se<-sd(dv)/sqrt(length(dv))
    } else {d<-mean(tapply(seq_along(sc),strata,function(i)mean(sc[i][grp[i]==1])-mean(sc[i][grp[i]==0])))
      se<-sqrt(var(sc[grp==1])/sum(grp==1)+var(sc[grp==0])/sum(grp==0))}
    res[[length(res)+1]]<-data.frame(study=st,module=mm,delta=d,se=se,n_hit=length(hit),n_bg=nrow(X))}}
out<-do.call(rbind,res); write.csv(out,"stage2/_orgdb_effects.csv",row.names=FALSE)
o<-out[out$module=="ISG_ISGF3",c("study","delta","se","n_bg")]
print(o,row.names=FALSE,digits=5)
