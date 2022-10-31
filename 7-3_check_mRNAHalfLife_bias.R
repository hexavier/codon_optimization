library(ggplot2)
library(ggpubr)

# Load data
TSprot = read.csv("results/TSall_PTR-UP-FCavg2-FCglob1.tsv", sep="\t")
halflife = read.delim("data/mRNA_HalfLife.tsv", sep="\t", row.names = 1)
cor(halflife, method = "spearman", use = "pairwise.complete.obs")

# For each tissue, check if there are differences in half-life
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t)&(TSprot$gene %in% rownames(halflife)),]
  dataset_temp = data.frame(row.names = seq(nrow(TStemp)*ncol(halflife)))
  dataset_temp$gene = rep(as.character(TStemp$gene),ncol(halflife))
  dataset_temp$halflife = as.numeric(unlist(halflife[as.character(TStemp$gene),]))
  dataset_temp$tissue = t
  dataset_temp$updown = rep(as.character(TStemp$updown),ncol(halflife))
  dataset_temp$HLcells = as.character(sapply(colnames(halflife), rep, nrow(TStemp)))
  dataset = rbind(dataset,dataset_temp)
}

# Plot
diffexp = compare_means(halflife~updown, data=dataset, group.by = c("tissue","HLcells"), method="wilcox.test")
ggplot(dataset, aes(x=tissue, y=halflife, fill=updown)) +
  facet_wrap(. ~ HLcells, ncol = 1, scales = "free") +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values=c("blue", "red")) +
  labs(x="Tissue", y = "mRNA Half-Life") +
  stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) +
  stat_compare_means(aes(label = ..p.adj..)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Get average over all 5 HL datasets
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t)&(TSprot$gene %in% rownames(halflife)),]
  dataset_temp = data.frame(row.names = seq(nrow(TStemp)))
  dataset_temp$gene = as.character(TStemp$gene)
  dataset_temp$halflife = rowMeans(halflife[as.character(TStemp$gene),],na.rm=T)
  dataset_temp$tissue = t
  dataset_temp$updown = as.character(TStemp$updown)
  dataset = rbind(dataset,dataset_temp)
}

# Plot
pdf("plots/mRNAHalfLife_avg_ALL.pdf",height=6,width =18)
diffexp = compare_means(halflife~updown, data=dataset, group.by = c("tissue"), method="wilcox.test")
print(ggplot(dataset, aes(x=tissue, y=halflife, fill=updown)) +
        geom_boxplot(position=position_dodge(1)) +
        labs(x="Tissue", y = "mRNA Half-Life") +
        scale_fill_manual(values=c("blue", "red")) +
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) +
        stat_compare_means(aes(label = ..p.adj..)) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
dev.off()