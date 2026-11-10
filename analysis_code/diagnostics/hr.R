suppressMessages({library(singscore);library(org.Mm.eg.db);library(org.Hs.eg.db);library(AnnotationDbi)})
CACHE<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
GS<-rbind(read.csv(file.path(CACHE,"genesets.csv"),stringsAsFactors=FALSE)[,c("set","genes")],
          read.csv(file.path(CACHE,"genesets_effector.csv"),stringsAsFactors=FALSE)[,c("set","genes")])
getg<-function(sn) unique(strsplit(GS$genes[GS$set==sn],";")[[1]])
res<-list()
score1<-function(dir,st,sfx){
  X<-as.matrix(read.csv(sprintf("%s/%s__expr.csv",dir,st),row.names=1,check.names=FALSE))
  md<-read.csv(sprintf("%s/%s__meta.csv",dir,st),colClasses="character")
  if(sfx=="_r"){rownames(X)<-toupper(rownames(X));sfx<-"_m"}
  rk<-rankGenes(X)
  for(m in c("ISG_ISGF3","EFF_TNFA_NFKB")){
    g<-getg(paste0(m,sfx)); if(sfx=="_m") g<-g
    hit<-intersect(g,rownames(X)); if(length(hit)<5) next
    s<-simpleScore(rk,upSet=hit,centerScore=TRUE)$TotalScore
    res[[length(res)+1]]<<-data.frame(study=st,module=m,
      ctrl_mean=mean(s[md$grp=="0"]), trt_mean=mean(s[md$grp=="1"]))}}
for(st in c("A_GSE112372_LPS3h","A_GSE112372_LPS8h")) score1("raw5",st,"_h")
score1("raw6","A_GSE110947_LPS4h","_h")
idx<-read.csv("stage2/_index.csv",stringsAsFactors=FALSE)
for(i in seq_len(nrow(idx))){
  sp<-idx$species[i]; sfx<-if(sp=="human")"_h" else if(sp=="rat")"_r" else "_m"
  try(score1("stage2",idx$study[i],sfx),silent=TRUE)}
out<-do.call(rbind,res); write.csv(out,"stage2/_headroom.csv",row.names=FALSE)
e<-out[out$module=="EFF_TNFA_NFKB",]
cat("EFF_TNFA_NFKB absolute singscore in CONTROL arms (centred score, range approx -0.5..+0.5)\n\n")
print(e[order(e$ctrl_mean),c("study","ctrl_mean","trt_mean")],row.names=FALSE,digits=3)
