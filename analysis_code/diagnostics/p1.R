suppressMessages({library(singscore);library(org.Mm.eg.db);library(org.Hs.eg.db);library(AnnotationDbi)})
C<-Sys.getenv("M1_CACHE", "cache_matrices")
G1<-read.csv(file.path(C,"genesets.csv"),stringsAsFactors=FALSE)
G2<-read.csv(file.path(C,"genesets_effector.csv"),stringsAsFactors=FALSE)
names(G2)[1]<-"set"; GS<-rbind(G1[,c("set","genes")],G2[,c("set","genes")])
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
cat("effector sets:",paste(G2$set,collapse=", "),"\n\n")
mk<-function(m,it,sc,sp){
  mode(m)<-"numeric"; m[!is.finite(m)]<-0
  L<-if(sc=="counts"){cs<-colSums(m);cs[cs==0]<-1;log2(t(t(m)/cs)*1e6+1)} else if(sc=="log2") m else log2(m+1)
  db<-if(sp=="human") org.Hs.eg.db else org.Mm.eg.db
  s<-if(it=="SYMBOL") rownames(L) else unname(suppressMessages(mapIds(db,keys=sub("\\..*$","",rownames(L)),keytype=it,column="SYMBOL",multiVals="first")))
  ok<-!is.na(s)&s!=""; L<-L[ok,,drop=FALSE]; rownames(L)<-s[ok]
  o<-order(rowMeans(L),decreasing=TRUE); L<-L[o,,drop=FALSE]; L[!duplicated(rownames(L)),,drop=FALSE] }
rd<-function(f,sep) read.csv(file.path(C,f),row.names=1,check.names=FALSE,sep=sep)
LD<-list(
 S1a_GSE274336=list("GSE274336_counts.csv.gz",",","ENSEMBL","counts","mouse"),
 S1b_GSE274019=list("GSE274019_20240603_counts.txt.gz","\t","ENTREZID","counts","mouse"),
 S2_GSE321887 =list("GSE321887_RNAseq_countMatrix.tsv.gz","\t","ENSEMBL","counts","human"),
 S4_GSE273834 =list("GSE273834_Salmon.norm.tpm.csv.gz",",","SYMBOL","tpm","mouse"),
 S5_GSE239377 =list("GSE239377_collagen_RNAseq_TPM.txt.gz","\t","ENSEMBL","tpm","mouse"),
 S7_GSE206771 =list("GSE206771_Counts.csv.gz",",","ENSEMBL","counts","mouse"),
 S8_GSE301422 =list("GSE301422_RawCounts.txt.gz","\t","ENSEMBL","counts","human"),
 S10_EMTAB7900=list("E-MTAB-7900_gene_counts_GRCr8.csv",",","SYMBOL","counts","mouse"),
 A_GSE274020_LPS6h=list("GSE274020_20240802_Raw_Counts.csv.gz",",","ENTREZID","counts","mouse"),
 A_GSE253423_LPS4h=list("GSE253423_raw_read_quant.table.txt.gz","\t","REFSEQ","counts","mouse"),
 A_GSE112372_LPS3h=list("GSE112372_raw_gene_counts.txt.gz","\t","SYMBOL","counts","human"))
MODS<-c(ISG="ISG_ISGF3",EFF_TNFA="EFF_TNFA_NFKB",EFF_INFL="EFF_INFLAMM",EFF_PANEL="EFF_FIELD_PANEL")
out<-list()
for(st in names(LD)){ L<-LD[[st]]
  M<-try(rd(L[[1]],L[[2]]),silent=TRUE); if(inherits(M,"try-error")){cat(st,"READ FAIL\n");next}
  m<-as.matrix(M); sp<-SP[SP$study==st,]
  # corrected contrasts
  if(st=="S4_GSE273834"){ sp<-sp[grepl("^WT",sp$sample),]; sp$stratum<-"a" }
  if(st=="S2_GSE321887"){ idx<-as.integer(sub(".*Sample_","",colnames(m))); pos<-((idx-1)%%3)+1; don<-((idx-1)%/%3)+1
    sel<-colnames(m)[pos %in% c(2,3)]
    sp<-data.frame(study=st,sample=sel,group=ifelse(((as.integer(sub(".*Sample_","",sel))-1)%%3)+1==3,"treated","control"),
                   stratum=paste0("D",((as.integer(sub(".*Sample_","",sel))-1)%/%3)+1),paired_id=NA,species="human",scale="counts") }
  cc<-intersect(sp$sample,colnames(m))
  if(length(cc)<nrow(sp)){cat(sprintf("%-18s sample mismatch %d/%d\n",st,length(cc),nrow(sp)));next}
  X<-mk(m[,cc,drop=FALSE],L[[3]],L[[4]],L[[5]]); R<-rankGenes(X)
  g<-ifelse(sp$group[match(cc,sp$sample)]=="treated",1,0); sv<-sp$stratum[match(cc,sp$sample)]
  sfx<-if(L[[5]]=="human") "_h" else "_m"
  row<-list(study=st,bg=nrow(X),n=length(cc))
  for(nm in names(MODS)){ key<-paste0(MODS[[nm]],sfx)
    gl<-unique(strsplit(GS$genes[GS$set==key],";")[[1]]); hit<-intersect(gl,rownames(X))
    if(length(hit)<5){row[[nm]]<-NA; next}
    sc<-simpleScore(R,upSet=hit,centerScore=TRUE)$TotalScore
    row[[nm]]<-if(length(unique(sv))>1) mean(sapply(unique(sv),function(u){k<-sv==u;mean(sc[k&g==1])-mean(sc[k&g==0])})) else mean(sc[g==1])-mean(sc[g==0]) }
  out[[st]]<-row
  cat(sprintf("%-18s bg=%6d n=%2d  ISG %+0.5f  TNFA %+0.5f  INFL %+0.5f  PANEL %+0.5f\n",
      st,row$bg,row$n,row$ISG,row$EFF_TNFA,row$EFF_INFL,row$EFF_PANEL)) }
saveRDS(out,"/tmp/p1.rds")
