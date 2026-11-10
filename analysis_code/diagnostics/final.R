suppressMessages({library(singscore);library(org.Mm.eg.db);library(org.Hs.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
G1<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
G2<-read.csv(file.path(C,"genesets_effector.csv"),stringsAsFactors=FALSE); names(G2)[1]<-"set"
GS<-rbind(G1[,c("set","genes")],G2[,c("set","genes")])
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
MODS<-c(ISG="ISG_ISGF3",TNFA="EFF_TNFA_NFKB",INFL="EFF_INFLAMM",PANEL="EFF_FIELD_PANEL")
mk<-function(m,it,sc,sp){ mode(m)<-"numeric"; m[!is.finite(m)]<-0
  L<-if(sc=="counts"){cs<-colSums(m);cs[cs==0]<-1;log2(t(t(m)/cs)*1e6+1)} else if(sc=="log2") m else log2(m+1)
  db<-if(sp=="human") org.Hs.eg.db else org.Mm.eg.db
  s<-if(it=="SYMBOL") rownames(L) else unname(suppressMessages(mapIds(db,keys=sub("\\..*$","",rownames(L)),keytype=it,column="SYMBOL",multiVals="first")))
  ok<-!is.na(s)&s!=""&s!="---"; L<-L[ok,,drop=FALSE]; rownames(L)<-s[ok]
  o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,,drop=FALSE]; L[!duplicated(rownames(L)),,drop=FALSE] }
rows<-list()
emit<-function(st,X,g,sv,pid,sp){ R<-rankGenes(X); sfx<-if(sp=="human")"_h" else "_m"
  for(nm in names(MODS)){ gl<-unique(strsplit(GS$genes[GS$set==paste0(MODS[[nm]],sfx)],";")[[1]])
    hit<-intersect(gl,rownames(X)); if(length(hit)<5) next
    sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
    if(!all(is.na(pid))&&length(unique(pid))>1&&!any(is.na(pid))){
      di<-sapply(unique(pid),function(p){k<-pid==p;sc[k&g==1]-sc[k&g==0]}); d<-mean(di); se<-sd(di)/sqrt(length(di))
    } else if(length(unique(sv))>1){
      di<-sapply(unique(sv),function(u){k<-sv==u;mean(sc[k&g==1])-mean(sc[k&g==0])}); d<-mean(di)
      se<-sqrt(var(sc[g==1])/sum(g==1)+var(sc[g==0])/sum(g==0))
    } else { d<-mean(sc[g==1])-mean(sc[g==0]); se<-sqrt(var(sc[g==1])/sum(g==1)+var(sc[g==0])/sum(g==0)) }
    rows[[length(rows)+1]]<<-data.frame(study=st,module=nm,species=sp,bg=nrow(X),n=length(g),
                                        nhit=length(hit),delta=d,se=se,stringsAsFactors=FALSE) } }
LD<-list(
 S1a_GSE274336=list("GSE274336_counts.csv.gz",",","ENSEMBL","counts","mouse"),
 S1b_GSE274019=list("GSE274019_20240603_counts.txt.gz","\t","ENTREZID","counts","mouse"),
 S2_GSE321887 =list("GSE321887_RNAseq_countMatrix.tsv.gz","\t","ENSEMBL","counts","human"),
 S4_GSE273834 =list("GSE273834_Salmon.norm.tpm.csv.gz",",","SYMBOL","tpm","mouse"),
 S5_GSE239377 =list("GSE239377_collagen_RNAseq_TPM.txt.gz","\t","ENSEMBL","tpm","mouse"),
 S7_GSE206771 =list("GSE206771_Counts.csv.gz",",","ENSEMBL","counts","mouse"),
 S8_GSE301422 =list("GSE301422_RawCounts.txt.gz","\t","ENSEMBL","counts","human"),
 S10_EMTAB7900=list("E-MTAB-7900_gene_counts_GRCr8.csv",",","SYMBOL","counts","mouse"),
 S11_GSE312997=list("GSE312997_sheet1_counts.csv",",","SYMBOL","counts","human"),
 A_GSE274020_LPS6h=list("GSE274020_20240802_Raw_Counts.csv.gz",",","ENTREZID","counts","mouse"),
 A_GSE253423_LPS4h=list("GSE253423_raw_read_quant.table.txt.gz","\t","REFSEQ","counts","mouse"),
 A_GSE112372_LPS3h=list("GSE112372_raw_gene_counts.txt.gz","\t","ENSEMBL","counts","human"))
for(st in names(LD)){ L<-LD[[st]]
  M<-try(read.csv(file.path(C,L[[1]]),row.names=1,check.names=FALSE,sep=L[[2]]),silent=TRUE)
  if(inherits(M,"try-error")) next
  m<-as.matrix(M); sp<-SP[SP$study==st,]
  if(st=="S4_GSE273834"){ sp<-sp[grepl("^WT",sp$sample),]; sp$stratum<-"a" }
  if(st=="S2_GSE321887"){ idx<-as.integer(sub(".*Sample_","",colnames(m))); pos<-((idx-1)%%3)+1
    sel<-colnames(m)[pos %in% c(2,3)]; si<-as.integer(sub(".*Sample_","",sel))
    sp<-data.frame(study=st,sample=sel,group=ifelse(((si-1)%%3)+1==3,"treated","control"),
                   stratum=paste0("D",((si-1)%/%3)+1),paired_id=paste0("D",((si-1)%/%3)+1),
                   species="human",scale="counts",stringsAsFactors=FALSE) }
  cc<-intersect(sp$sample,colnames(m)); if(length(cc)<nrow(sp)) next
  i<-match(cc,sp$sample)
  emit(st,mk(m[,cc,drop=FALSE],L[[3]],L[[4]],L[[5]]),ifelse(sp$group[i]=="treated",1,0),sp$stratum[i],sp$paired_id[i],L[[5]]) }
# S6 from tar
fs<-sort(list.files("/tmp/g291",full.names=TRUE)); gsm<-sub("_.*","",basename(fs))
ms<-lapply(fs,function(f) read.csv(f,sep="\t",header=TRUE,row.names=1,check.names=FALSE)[,1])
m6<-do.call(cbind,ms); rownames(m6)<-rownames(read.csv(fs[1],sep="\t",header=TRUE,row.names=1,check.names=FALSE)); colnames(m6)<-gsm
sp<-SP[SP$study=="S6_GSE291556",]; cc<-intersect(sp$sample,colnames(m6)); i<-match(cc,sp$sample)
emit("S6_GSE291556",mk(m6[,cc,drop=FALSE],"SYMBOL","norm","mouse"),ifelse(sp$group[i]=="treated",1,0),rep("a",length(cc)),rep(NA,length(cc)),"mouse")
# S3 array
ln<-readLines(gzfile(file.path(C,"GSE36878_series_matrix.txt.gz")))
a<-grep("^!series_matrix_table_begin",ln); b<-grep("^!series_matrix_table_end",ln)
tb<-read.csv(text=paste(ln[(a+1):(b-1)],collapse="\n"),sep="\t",row.names=1,check.names=FALSE)
an<-read.csv(file.path(C,"GPL11180_annot.tsv"),sep="\t",check.names=FALSE,quote="")
map<-setNames(as.character(an$symbol),as.character(an$ID))
mm<-as.matrix(tb); mode(mm)<-"numeric"; mm[!is.finite(mm)]<-0
Lx<-if(max(mm,na.rm=TRUE)>50) log2(mm+1) else mm
s<-unname(map[rownames(Lx)]); ok<-!is.na(s)&s!=""&s!="---"; Lx<-Lx[ok,,drop=FALSE]; rownames(Lx)<-s[ok]
o<-order(rowMeans(Lx),decreasing=TRUE); Lx<-Lx[o,]; Lx<-Lx[!duplicated(rownames(Lx)),]
sp<-SP[SP$study=="S3_GSE36878",]; cc<-intersect(sp$sample,colnames(Lx)); i<-match(cc,sp$sample)
emit("S3_GSE36878",Lx[,cc,drop=FALSE],ifelse(sp$group[i]=="treated",1,0),sp$stratum[i],rep(NA,length(cc)),"mouse")
D<-do.call(rbind,rows)
write.csv(D,"/tmp/unfiltered_effects.csv",row.names=FALSE)
cat(sprintf("wrote %d rows, %d studies\n",nrow(D),length(unique(D$study))))
print(D[D$module=="ISG",c("study","bg","n","delta","se")],row.names=FALSE)
