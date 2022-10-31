### PTR ###
## Load GTEX data
gtex = read.csv("data/PTR_GTEx_Jiang2020.csv",row.names = 1)
mapping = read.csv("data/GTEx_samples.csv")

# Map to tissue and get median of GTEx
tissues = as.character(sapply(colnames(gtex), function(x) mapping[mapping$Sample.ID==x,"Tissue"]))
tis_med = sapply(unique(tissues), function(x) apply(gtex[,tissues %in% x, drop=F],1,median, na.rm=T))

## Load Eraslan PTR
Eraslan = read.csv("data/PTR_Eraslan2019.tsv", sep="\t", row.names = 2)[,4:32]

## Correlation PTR between replicates of GTEx
within_tis = sapply(unique(tissues), function(x) cor(gtex[,tissues %in% x, drop=F],method="spearman",use="pairwise"))
within_tis = unlist(lapply(within_tis,function(x) mean(x[upper.tri(x)])))

## Correlation PTR between non-matched tissues
ENScons = rownames(gtex)[rownames(gtex) %in% rownames(Eraslan)]
matched = colnames(Eraslan)[colnames(Eraslan) %in% unique(tissues)]
nonmatching_tis_gtex = sapply(matched, function(x) mean(cor(tis_med[ENScons,x],tis_med[ENScons,!(colnames(tis_med) %in% x)],method="spearman",use="pairwise")))
nonmatching_tis_eraslan = sapply(matched, function(x) mean(cor(Eraslan[ENScons,x],Eraslan[ENScons,!(colnames(Eraslan) %in% x)],method="spearman",use="pairwise")))

## Correlation PTR along matching tissues (for each gene)
matching_tis = sapply(ENScons, function(x) if(sum(colSums(is.na(rbind(tis_med[x,matched],Eraslan[x,matched])))<1)>5){
  cor(as.numeric(tis_med[x,matched]),as.numeric(Eraslan[x,matched]),method="spearman",use="pairwise")}else{NA})
hist(matching_tis,breaks=20)
# Permutation test
outtest = c()
for (n in seq(1000)){
  rand_matched = sample(matched)
  random_tis = sapply(ENScons, function(x) if(sum(colSums(is.na(rbind(tis_med[x,matched],Eraslan[x,matched])))<1)>5){
    cor(as.numeric(tis_med[x,matched]),as.numeric(Eraslan[x,rand_matched]),method="spearman",use="pairwise")}else{NA})
  outtest = c(outtest,median(matching_tis - random_tis, na.rm=T)>0)
}
sprintf("Permutation test p-val: %f", sum(outtest==F)/length(outtest))

### mRNA ###
mrnaeraslan = read.csv("data/mRNAseq_FPKM_Eraslan2019.csv",row.names = 2)[,5:33]

mrnagtex = read.csv("data/mRNAseq_GTEx_Jiang2020.csv",row.names = 1,skip = 1)
# Remove TPM<10 (same as Eraslan)
mrnagtex[mrnagtex<log2(10)] = NA; mrnagtex = mrnagtex[rowSums(!is.na(mrnagtex))>0,]
mrna_samples = sapply(colnames(mrnagtex),function(x) paste0(unlist(strsplit(x,"\\."))[1:3],collapse="."))
colnames(mrnagtex) = mrna_samples
tissues = as.character(sapply(colnames(mrnagtex), function(x) mapping[mapping$Sample.ID==x,"Tissue"]))
mrna_med = sapply(unique(tissues), function(x) apply(mrnagtex[,tissues %in% x, drop=F],1,median, na.rm=T))
write.csv(mrna_med/log2(10),"data/mRNAseq_medians_GTEx_Jiang2020.csv")

ENScons_mrna = rownames(mrnaeraslan)[rownames(mrnaeraslan) %in% rownames(mrna_med)]
mrna_cor = sapply(ENScons_mrna, function(x) if(sum(colSums(is.na(rbind(mrna_med[x,matched],mrnaeraslan[x,matched])))<1)>5){
  cor(as.numeric(mrna_med[x,matched]),as.numeric(mrnaeraslan[x,matched]), method="spearman",use="pairwise")}else{NA})
hist(mrna_cor,breaks=20)
# Permutation test
outtest = c()
for (n in seq(1000)){
  rand_matched = sample(matched)
  random_cor = sapply(ENScons_mrna, function(x) if(sum(colSums(is.na(rbind(mrna_med[x,matched],mrnaeraslan[x,matched])))<1)>5){
    cor(as.numeric(mrna_med[x,matched]),as.numeric(mrnaeraslan[x,rand_matched]), method="spearman",use="pairwise")}else{NA})
  outtest = c(outtest,median(mrna_cor - random_cor, na.rm=T)>0)
}
sprintf("Permutation test p-val: %f", sum(outtest==F)/length(outtest))

### Prot ###
proteraslan = read.csv("data/proteome_iBAQ_Eraslan2019.csv", row.names = 2)[,5:33]

protgtex = read.csv("data/proteome_GTEx_Jiang2020.csv",row.names = 1, skip = 1)
# Remove 0s because of imputation
protgtex[protgtex==0] = NA; protgtex = protgtex[rowSums(!is.na(protgtex))>0,]
prot_samples = sapply(colnames(protgtex),function(x) paste0(unlist(strsplit(x,"\\."))[1:3],collapse="."))
colnames(protgtex) = prot_samples
tissues = as.character(sapply(colnames(protgtex), function(x) mapping[mapping$Sample.ID==x,"Tissue"]))
prot_med = sapply(unique(tissues), function(x) apply(protgtex[,tissues %in% x, drop=F],1,median, na.rm=T))
write.csv(prot_med/log2(10),"data/proteome_medians_GTEx_Jiang2020.csv")

ENScons_prot = rownames(proteraslan)[rownames(proteraslan) %in% rownames(prot_med)]
prot_cor = sapply(ENScons_prot, function(x) if(sum(colSums(is.na(rbind(prot_med[x,matched],proteraslan[x,matched])))<1)>5){
  cor(as.numeric(prot_med[x,matched]),as.numeric(proteraslan[x,matched]),method="spearman",use="pairwise")}else{NA})
hist(prot_cor,breaks=20)
# Permutation test
outtest = c()
for (n in seq(1000)){
  rand_matched = sample(matched)
  random_cor = sapply(ENScons_prot, function(x) if(sum(colSums(is.na(rbind(prot_med[x,matched],proteraslan[x,matched])))<1)>5){
    cor(as.numeric(prot_med[x,matched]),as.numeric(proteraslan[x,rand_matched]),method="spearman",use="pairwise")}else{NA})
  outtest = c(outtest,median(prot_cor - random_cor, na.rm=T)>0)
}
sprintf("Permutation test p-val: %f", sum(outtest==F)/length(outtest))


### COMPARE VARIABILITY OF PROT VS mRNA ###
sdgtex = mean(apply(prot_med[ENScons_prot,matched],1,sd,na.rm=F),na.rm=T)/mean(apply(mrna_med[ENScons_mrna,matched],1,sd,na.rm=F),na.rm=T)
sderaslan = mean(apply(proteraslan[ENScons_prot,matched],1,sd,na.rm=F),na.rm=T)/mean(apply(mrnaeraslan[ENScons_mrna,matched],1,sd,na.rm=F),na.rm=T)
sderaslan/sdgtex

### Prot vs mRNA ###
ENScons_all = ENScons_mrna[ENScons_mrna %in% ENScons_prot]
eraslan_cor = sapply(matched, function(x) cor(mrnaeraslan[ENScons_all,x],proteraslan[ENScons_all,x],method="spearman",use="pairwise"))
