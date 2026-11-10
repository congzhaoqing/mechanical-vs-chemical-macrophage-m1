suppressMessages({library(singscore)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
G1<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
G2<-read.csv(file.path(C,"genesets_effector.csv"),stringsAsFactors=FALSE); names(G2)[1]<-"set"
GS<-rbind(G1[,c("set","genes")],G2[,c("set","genes")])
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
ln<-readLines(gzfile(file.path(C,"GSE36878_series_matrix.txt.gz")))
st<-grep("^!series_matrix_table_begin",ln); en<-grep("^!series_matrix_table_end",ln)
tb<-read.csv(text=paste(ln[(st+1):(en-1)],collapse="\n"),sep="\t",row.names=1,check.names=FALSE)
plat<-sub('.*\t','',grep("^!Series_platform_id",ln,value=TRUE)[1]); plat<-gsub('"','',plat)
cat("platform:",plat," matrix:",nrow(tb),"x",ncol(tb),"\n")
an<-read.csv(file.path(C,paste0(plat,"_annot.tsv")),sep="\t",check.names=FALSE,quote="")
cat("annot cols:",paste(head(names(an),4),collapse=","),"\n")
idc<-names(an)[1]; symc<-grep("symbol",names(an),ignore.case=TRUE,value=TRUE)[1]
map<-setNames(as.character(an[[symc]]),as.character(an[[idc]]))
s<-unname(map[rownames(tb)])
m<-as.matrix(tb); mode(m)<-"numeric"; m[!is.finite(m)]<-0
L<-if(max(m,na.rm=TRUE)>50) log2(m+1) else m
ok<-!is.na(s)&s!=""&s!="---"; L<-L[ok,,drop=FALSE]; rownames(L)<-s[ok]
o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,]; L<-L[!duplicated(rownames(L)),]
sp<-SP[SP$study=="S3_GSE36878",]; cc<-intersect(sp$sample,colnames(L))
cat(sprintf("S3 matched %d/%d samples\n",length(cc),nrow(sp)))
if(length(cc)==nrow(sp)){ X<-L[,cc,drop=FALSE]; R<-rankGenes(X)
  g<-ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0); sv<-sp$stratum[match(cc,sp$sample)]
  v<-c(); for(nm in c("ISG_ISGF3_m","EFF_TNFA_NFKB_m","EFF_INFLAMM_m","EFF_FIELD_PANEL_m")){
    gl<-unique(strsplit(GS$genes[GS$set==nm],";")[[1]]); hit<-intersect(gl,rownames(X))
    sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
    v[nm]<-if(length(unique(sv))>1) mean(sapply(unique(sv),function(u){k<-sv==u;mean(sc[k&g==1])-mean(sc[k&g==0])})) else mean(sc[g==1])-mean(sc[g==0]) }
  cat(sprintf("S3_GSE36878        bg=%6d n=%2d  ISG %+0.5f  TNFA %+0.5f  INFL %+0.5f  PANEL %+0.5f   (published ISG -0.000223)\n",
      nrow(X),length(cc),v[1],v[2],v[3],v[4]))
  saveRDS(v,"/tmp/p3.rds") }
