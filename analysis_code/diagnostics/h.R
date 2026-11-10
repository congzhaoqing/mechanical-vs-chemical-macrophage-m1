suppressMessages(library(singscore))
CACHE<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS<-rbind(read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)[,c("set","genes")],
          read.csv(file.path(CACHE,"genesets_effector.csv"),stringsAsFactors=FALSE)[,c("set","genes")])
getg<-function(sn) unique(strsplit(GS$genes[GS$set==sn],";")[[1]])
MODS<-c("ISG_ISGF3","ISG_IFNA","ISG_IFNG","YAP","EFF_TNFA_NFKB","EFF_INFLAMM","EFF_FIELD_PANEL")
idx<-read.csv("raw4/_index.csv",stringsAsFactors=FALSE); res<-list()
for(i in seq_len(nrow(idx))){ st<-idx$study[i]
  X<-as.matrix(read.csv(sprintf("raw4/%s__expr.csv",st),row.names=1,check.names=FALSE))
  md<-read.csv(sprintf("raw4/%s__meta.csv",st),colClasses="character")
  grp<-as.integer(md$grp); paired<-md$paired; rk<-rankGenes(X)
  for(m in MODS){ hit<-intersect(getg(paste0(m,"_h")),rownames(X)); if(length(hit)<5) next
    s<-simpleScore(rk,upSet=hit,centerScore=TRUE)$TotalScore
    if(any(paired!="")){dv<-tapply(seq_along(s),paired,function(i)mean(s[i][grp[i]==1])-mean(s[i][grp[i]==0]))
      d<-mean(dv,na.rm=TRUE);se<-sd(dv,na.rm=TRUE)/sqrt(sum(!is.na(dv)))
    } else {d<-mean(s[grp==1])-mean(s[grp==0])
      se<-sqrt(var(s[grp==1])/sum(grp==1)+var(s[grp==0])/sum(grp==0))}
    res[[length(res)+1]]<-data.frame(study=st,module=m,delta=d,se=se,n_hit=length(hit))}}
out<-do.call(rbind,res); write.csv(out,"raw4/_effects.csv",row.names=FALSE)
w<-reshape(out[,c("study","module","delta")],idvar="study",timevar="module",direction="wide")
names(w)<-sub("delta.","",names(w)); print(w,row.names=FALSE,digits=3)
