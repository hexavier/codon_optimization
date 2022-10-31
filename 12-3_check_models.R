library(ggplot2)

# Load proteome data
proteraslan = read.csv("data/proteome_iBAQ_Eraslan2019.csv", row.names = 1)
protmann_raw = read.csv("data/proteome_iBAQ_Mann2012.csv") # from http://maxqb.biochem.mpg.de/mxdb/
# Unique IDs
uniqid = unique(as.character(protmann_raw$Leading.Gene.Name)); uniqid = uniqid[uniqid!=""]
protmann = t(sapply(uniqid, function(x) apply(protmann_raw[protmann_raw$Leading.Gene.Name==x,c("Intensity.A549_1","Intensity.A549_2","Intensity.A549_3","Intensity.HEK293_1","Intensity.HEK293_2","Intensity.HEK293_3"),drop=F],2,max,na.rm=T)))

# Map to ENSG id
rownames(protmann) = sapply(rownames(protmann),function(x) if(x %in% rownames(proteraslan)){proteraslan[x,"EnsemblGeneID"]}else{x})
rownames(proteraslan) = proteraslan$EnsemblGeneID
proteraslan = proteraslan[,5:33]

# Remove 0 to NA and log
protmann[protmann==0] <- NA
protmann = log10(protmann)
# Take median
protmann_med = sapply(c("A549","HEK293"),function(x) apply(protmann[,grepl(x,colnames(protmann))],1,median,na.rm=T))
protmann_med = protmann_med[rowSums(!is.na(protmann_med))==2,]

# Correlate cell lines and compute differences
commonprot = rownames(proteraslan)[rownames(proteraslan) %in% rownames(protmann_med)]
cortab = as.matrix(cor(cbind(proteraslan[commonprot,],protmann_med[commonprot,]),method = "spearman", use="pairwise"))[1:29,30:31]
normcortab = as.data.frame(apply(cortab,2,scale)); rownames(normcortab) = rownames(cortab)
normcortab$diff = normcortab$A549 - normcortab$HEK293
normcortab$tissue = factor(rownames(normcortab),levels=rownames(normcortab)[order(normcortab$diff,decreasing = T)])

# Plot
ggplot(normcortab, aes(x=tissue, y=diff)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y="Cor_A549 - Cor_HEK293")
ggsave("plots/cell_lines_correlations.pdf",width = 5.5, height = 2.5)

cortab = as.data.frame(cortab); cortab$lab = rownames(cortab)
ggplot(cortab, aes(x=A549, y=HEK293, label=lab)) +
  geom_point() +
  geom_text(data=subset(cortab, lab %in% c("Kidney", "Lung"))) +
  theme_classic() +
  labs(x="Cor_A549", y="Cor_HEK293")
ggsave("plots/cell_lines_correlations_scatter.pdf",width = 2, height = 2)

### Compare similarity of codon usage ###
codus = read.csv("data/human_CU_ENS_09062020.tsv",sep="\t")[,5:68]

# Make conformable matrices
consENS = commonprot[commonprot %in% rownames(codus)]
protmerged = cbind(proteraslan[consENS,],protmann_med[consENS,])
protmerged[is.na(protmerged)] <- 0
weightCU = t(protmerged) %*% as.matrix(codus[consENS,])
normCU = apply(weightCU,1,function(x) x/sum(x))

# Correlation
CUcortab = as.matrix(cor(normCU,method = "pearson", use="pairwise"))[1:29,30:31]
CUnormcortab = as.data.frame(apply(CUcortab,2,scale)); rownames(CUnormcortab) = rownames(CUcortab)
CUnormcortab$diff = CUnormcortab$A549 - CUnormcortab$HEK293
CUnormcortab$tissue = factor(rownames(CUnormcortab),levels=rownames(CUnormcortab)[order(CUnormcortab$diff,decreasing = T)])

# Plot
ggplot(CUnormcortab, aes(x=tissue, y=diff)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y="Cor_A549 - Cor_HEK293")
ggsave("plots/cell_lines_codus_correlations.pdf",width = 5.5, height = 2.5)

CUcortab = as.data.frame(CUcortab); CUcortab$lab = rownames(CUcortab)
ggplot(CUcortab, aes(x=A549, y=HEK293, label=lab)) +
  geom_point() +
  geom_text(data=subset(CUcortab, lab %in% c("Kidney", "Lung"))) +
  theme_classic() +
  labs(x="Cor_A549", y="Cor_HEK293")
ggsave("plots/cell_lines_codus_correlations_scatter.pdf",width = 2, height = 2)
