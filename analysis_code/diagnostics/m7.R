suppressMessages(library(singscore))
CACHE<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-rbind(read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)[,c("set","genes")],
          read.csv(file.path(CACHE,"genesets_effector.csv"),stringsAsFactors=FALSE)[,c("set","genes")])
getg<-function(sn) unique(strsplit(GS$genes[GS$set==sn],";")[[1]])
MODS<-c("ISG_ISGF3","ISG_IFNA","ISG_IFNG","YAP","EFF_TNFA_NFKB","EFF_INFLAMM","EFF_FIELD_PANEL")
idx<-read.csv("stage7/_index.csv",stringsAsFactors=FALSE); res<-list()
for(i in seq_len(nrow(idx))){
  st<-idx$study[i]; sfx<-if(idx$species[i]=="human")"_h" else "_m"
  X<-as.matrix(read.csv(sprintf("stage7/%s__expr.csv",st),row.names=1,check.names=FALSE))
  md<-read.csv(sprintf("stage7/%s__meta.csv",st),colClasses="character")
  grp<-as.integer(md$grp); paired<-ifelse(is.na(md$paired),"",md$paired); rk<-rankGenes(X)
  for(m in MODS){ hit<-intersect(getg(paste0(m,sfx)),rownames(X)); if(length(hit)<5) next
    s<-simpleScore(rk,upSet=hit,centerScore=TRUE)$TotalScore
    if(any(paired!="")){dv<-tapply(seq_along(s),paired,function(i)mean(s[i][grp[i]==1])-mean(s[i][grp[i]==0]))
      d<-mean(dv,na.rm=TRUE); se<-sd(dv,na.rm=TRUE)/sqrt(sum(!is.na(dv)))
    } else {d<-mean(s[grp==1])-mean(s[grp==0])
      se<-if(sum(grp==1)>1 && sum(grp==0)>1) sqrt(var(s[grp==1])/sum(grp==1)+var(s[grp==0])/sum(grp==0)) else NA}
    res[[length(res)+1]]<-data.frame(study=st,module=m,delta=d,se=se,n_hit=length(hit),n_bg=nrow(X))}}
out<-do.call(rbind,res); write.csv(out,"stage7/_effects.csv",row.names=FALSE)
cat("rows",nrow(out),"| contrasts",length(unique(out$study)),"\n")
