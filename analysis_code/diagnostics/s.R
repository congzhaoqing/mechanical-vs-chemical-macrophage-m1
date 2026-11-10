C<-"/Users/zhaoqingcong/Documents/My Mac/Documents/Stony Brook/Code/Claude/bioinfo project/执行v5/cache_matrices"
SP<-read.csv(file.path(C,"_contrast_spec.csv"),stringsAsFactors=FALSE)
for(f in c("GSE321887_RNAseq_countMatrix.tsv.gz","GSE273834_Salmon.norm.tpm.csv.gz","GSE239377_collagen_RNAseq_TPM.txt.gz")){
  sep<-if(grepl("tsv|txt",f)) "\t" else ","
  M<-read.csv(file.path(C,f),row.names=1,check.names=FALSE,sep=sep,nrows=2)
  cat("\n===",f,"\n  ALL columns (",ncol(M),"):\n"); print(colnames(M))
}
cat("\n=== spec groupings used ===\n")
for(st in c("S2_GSE321887","S4_GSE273834","S5_GSE239377")){
  s<-SP[SP$study==st,]; cat("\n",st,"\n"); print(data.frame(sample=s$sample,group=s$group,stratum=s$stratum))
}
