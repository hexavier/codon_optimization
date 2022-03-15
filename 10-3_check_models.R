library(ggplot2)

# Load proteome data
proteraslan = read.csv("data/proteome_iBAQ_Eraslan2019.csv", row.names = 1)[,5:33]
protmann_raw = read.csv("data/proteome_iBAQ_Mann2012.csv") # from http://maxqb.biochem.mpg.de/mxdb/

# Unique IDs
uniqid = unique(as.character(protmann_raw$Leading.Gene.Name)); uniqid = uniqid[uniqid!=""]
protmann = t(sapply(uniqid, function(x) apply(protmann_raw[protmann_raw$Leading.Gene.Name==x,c("Intensity.A549_1","Intensity.A549_2","Intensity.A549_3","Intensity.HEK293_1","Intensity.HEK293_2","Intensity.HEK293_3"),drop=F],2,max,na.rm=T)))
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(y="Cor_A549 - Cor_HEK293")
ggsave("plots/cell_lines_correlations.pdf",width = 5, height = 3)
          