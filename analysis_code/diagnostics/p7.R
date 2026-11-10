suppressMessages(library(singscore))
CACHE<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS<-rbind(read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)[,c("set","genes")],
          read.csv(file.path(CACHE,"genesets_effector.csv"),stringsAsFactors=FALSE)[,c("set","genes")])
getg<-function(sn) unique(strsplit(GS$genes[GS$set==sn],";")[[1]])
MODS<-c("ISG_ISGF3","EFF_TNFA_NFKB")
dof<-function(s,g,p){ if(any(p!="")) mean(tapply(seq_along(s),p,function(i)mean(s[i][g[i]==1])-mean(s[i][g[i]==0])),na.rm=TRUE)
                      else mean(s[g==1])-mean(s[g==0]) }
idx<-read.csv("stage7/_index.csv",stringsAsFactors=FALSE); res<-list()
set.seed(42)
for(i in seq_len(nrow(idx))){
  st<-idx$study[i]; sfx<-if(idx$species[i]=="human")"_h" else "_m"
  X<-as.matrix(read.csv(sprintf("stage7/%s__expr.csv",st),row.names=1,check.names=FALSE))
  md<-read.csv(sprintf("stage7/%s__meta.csv",st),colClasses="character")
  g<-as.integer(md$grp); p<-ifelse(is.na(md$paired),"",md$paired); rk<-rankGenes(X); n<-length(g)
  if(any(p!="")){ up<-unique(p); perms<-as.matrix(expand.grid(rep(list(c(0,1)),length(up))))
    P<-t(apply(perms,1,function(fl){ v<-g; for(j in seq_along(up)) if(fl[j]==1) v[p==up[j]]<-1-v[p==up[j]]; v }))
  } else { cmb<-combn(n,sum(g==1),simplify=FALSE)
    P<-t(sapply(cmb,function(k){v<-rep(0,n);v[k]<-1;v})) }
  if(nrow(P)>20000) P<-P[sample(nrow(P),20000),]
  for(m in MODS){ hit<-intersect(getg(paste0(m,sfx)),rownames(X)); if(length(hit)<5) next
    s<-simpleScore(rk,upSet=hit,centerScore=TRUE)$TotalScore
    d<-dof(s,g,p); pd<-apply(P,1,function(gg) dof(s,gg,p))
    res[[length(res)+1]]<-data.frame(study=st,module=m,delta=d,
      p_perm=mean(abs(pd)>=abs(d)-1e-12), perm_space=nrow(P))}}
out<-do.call(rbind,res); write.csv(out,"stage7/_perm.csv",row.names=FALSE)
cat("done",nrow(out),"rows\n")
