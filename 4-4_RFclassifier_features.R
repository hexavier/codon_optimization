library(ggplot2)
library(ggpubr)
library(gplots)
library(ggrepel)

# Load data
codus_UPDOWN = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)[,3:62]
rcu = read.csv("results/mean_all_randomforest_RCU_to_UPDOWNproteins.csv", row.names = 1)[,3:62]
codusnormratios = read.csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
RCUratios = read.csv("results/all_RCU_ratios_UPDOWNproteins.csv", row.names = 1)

# Define structure
codons = colnames(codus_UPDOWN); codons = codons[-4]
dataset = c()
for (t in rownames(codus_UPDOWN)){
  dataset_temp = data.frame(row.names = 1:length(codons))
  dataset_temp$codus_UPDOWN = as.numeric(codus_UPDOWN[t,codons])
  dataset_temp$rcu = as.numeric(rcu[t,codons])
  dataset_temp$codusnormratios = as.numeric(codusnormratios[t,codons])
  dataset_temp$RCUratios = as.numeric(RCUratios[t,codons])
  dataset_temp$tissue = t
  dataset_temp$codon = codons
  # Merge
  dataset = rbind(dataset,dataset_temp)
}

pdf("plots/ALL_codusnorm_vs_othermetrics.pdf",height=15,width =20)

print(ggplot(dataset, aes(x=codus_UPDOWN, y=abs(codusnormratios))) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=codus_UPDOWN, y=abs(RCUratios))) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=rcu, y=abs(RCUratios))) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=codus_UPDOWN, y=rcu)) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=codusnormratios, y=RCUratios)) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
dev.off()

# Idenitfy differences tissues
outweights <- heatmap.2(t(codus_UPDOWN),symm=F, margins = c(6,5), na.color = "grey", trace="none")
cor(t(codus_UPDOWN), method = "spearman", use = "pairwise")
outratios <- heatmap.2(t(codusnormratios),col=bluered, symm=F, margins = c(6,5), na.color = "grey", trace="none")
cor(t(codusnormratios), method = "spearman", use = "pairwise")

# Barplot of AUC
tissueorder = rownames(codus_UPDOWN)[outweights$colInd]
aucs = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)[tissueorder,1]
barplot(aucs,ylim=c(0,1)); abline(h=0.5)

tissueorder = rownames(codusnormratios)[outratios$colInd]
aucs = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)[tissueorder,1]
barplot(aucs,ylim=c(0,1)); abline(h=0.5)

#### CHECK ASSOCIATION WITH Ki67 ####
# Get clusters
hc <- hclust(dist(codusnormratios))
clusters = cutree(hc, k=2)
at = rownames(codusnormratios)[clusters==1]
gc = rownames(codusnormratios)[clusters==2]
# Load expression data
protgtex = read.csv("data/proteome_medians_GTEx_Jiang2020.csv",row.names = 1)
# Build dataset
expdata = as.data.frame(t(protgtex["ENSG00000148773",,drop=F]))
expdata$cluster = sapply(rownames(expdata), function(x) if(x %in% at){"AT Cluster"}else{"GC Cluster"})
expdata$tissue = rownames(expdata)
ggplot(expdata, aes(x=cluster, y=ENSG00000148773, color=cluster, label=tissue)) +
  geom_point(position=position_jitter(w = 0.1, h = 0)) +
  labs(y = "Ki67 Protein Abundance") +
  scale_color_manual(values=c("red", "blue")) +
  stat_compare_means(method = "t.test") +
  geom_text_repel(size=3) +
  theme_classic()
ggsave("plots/association_ki67.pdf",width = 4, height = 4)

