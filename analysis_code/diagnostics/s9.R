suppressMessages(library(singscore))
CACHE<-Sys.getenv("M1_CACHE", "cache_matrices")
GS<-rbind(read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)[,c("set","genes")],
          read.csv(file.path(CACHE,"genesets_effector.csv"),stringsAsFactors=FALSE)[,c("set","genes")])
getg<-function(sn) unique(strsplit(GS$genes[GS$set==sn],";")[[1]])
MODS<-c("ISG_ISGF3","EFF_TNFA_NFKB","EFF_INFLAMM","EFF_FIELD_PANEL","YAP")
idx<-read.csv("stage8/_index.csv",stringsAsFactors=FALSE); res<-list(); set.seed(42)
for(i in seq_len(nrow(idx))){
  st<-idx$study[i]
  X<-as.matrix(read.csv(sprintf("stage8/%s__expr.csv",st),row.names=1,check.names=FALSE))
  md<-read.csv(sprintf("stage8/%s__meta.csv",st),colClasses="character")
  g<-as.integer(md$grp); sr<-md$stratum; n<-length(g); rk<-rankGenes(X)
  dof<-function(s,gg) mean(tapply(seq_along(s),sr,function(i)mean(s[i][gg[i]==1])-mean(s[i][gg[i]==0])),na.rm=TRUE)
  if(length(unique(sr))>1){
    P<-as.matrix(expand.grid(rep(list(c(0,1)),length(unique(sr)))))
    PM<-t(apply(P,1,function(fl){v<-g;for(j in seq_along(unique(sr))){k<-unique(sr)[j]; if(fl[j]==1) v[sr==k]<-1-v[sr==k]};v}))
  } else { cmb<-combn(n,sum(g==1),simplify=FALSE); PM<-t(sapply(cmb,function(k){v<-rep(0,n);v[k]<-1;v})) }
  bg<-rownames(X)
  for(m in MODS){ hit<-intersect(getg(paste0(m,"_m")),bg); if(length(hit)<5) next
    s<-simpleScore(rk,upSet=hit,centerScore=TRUE)$TotalScore; d<-dof(s,g)
    pd<-apply(PM,1,function(gg) dof(s,gg)); pp<-mean(abs(pd)>=abs(d)-1e-12)
    rnd<-replicate(500,dof(simpleScore(rk,upSet=sample(bg,length(hit)),centerScore=TRUE)$TotalScore,g))
    res[[length(res)+1]]<-data.frame(study=st,module=m,delta=d,p_perm=pp,perm_space=nrow(PM),
      z_rand=(d-mean(rnd))/sd(rnd), n_hit=length(hit), n_bg=length(bg))}}
out<-do.call(rbind,res); write.csv(out,"stage8/_effects.csv",row.names=FALSE)
cat("rows",nrow(out),"\n")
