## Load data
mrna = read.csv("data/mRNAseq_GTEx_Jiang2020.csv",row.names = 1,skip = 1)
prot = read.csv("data/proteome_GTEx_Jiang2020.csv",row.names = 1, skip = 1)
# Remove 0s because of imputation and remove TPM<10 (same as Eraslan)
prot[prot==0] = NA; mrna[mrna<log2(10)] = NA
codus = read.csv("data/human_CU_ENS_09062020.tsv", row.names = 1, sep = "\t")
mapping = read.csv("data/GTEx_samples.csv")

## Make median ctt of protein quantifications (same as Eraslan)
overall = median(as.matrix(prot), na.rm=T)
prot = apply(prot,2,function(x) x - (median(x,na.rm=T)-overall))

## Compute p/m ratios per tissue
# Consensus genes
ENScons = rownames(codus)[rownames(codus) %in% rownames(mrna)[rownames(mrna) %in% rownames(prot)]]
# Consensus samples
mrna_samples = sapply(colnames(mrna),function(x) paste0(unlist(strsplit(x,"\\."))[1:3],collapse="."))
colnames(mrna) = mrna_samples
prot_samples = sapply(colnames(prot),function(x) paste0(unlist(strsplit(x,"\\."))[1:3],collapse="."))
colnames(prot) = prot_samples
GTEXcons = prot_samples[prot_samples %in% mrna_samples]

# Ratios
# As data is log-normalized, the ratio corresponds to a subtraction. Change base 2 to base 10
ratios = prot[ENScons,GTEXcons] - mrna[ENScons,GTEXcons]
ratios = ratios/log2(10)
write.csv(ratios,"data/PTR_GTEx_Jiang2020.csv")

# Map to tissue and get median
tissues = as.character(sapply(GTEXcons, function(x) mapping[mapping$Sample.ID==x,"Tissue"]))
tis_med = sapply(unique(tissues), function(x) apply(ratios[,tissues %in% x, drop=F],1,median, na.rm=T))

## Save
write.csv(tis_med,"data/PTR_medians_GTEx_Jiang2020.csv")
