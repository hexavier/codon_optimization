library(ggplot2)
library(ggpubr)

# Load data
TSprot = read.csv("results/TSall_PTR-UP-FCavg2-FCglob1.tsv", sep="\t")
Kuster = read.csv("data/KusterProteinHalfLife.csv")
Mathieson = read.csv("data/MathiesonProteinHalfLife.csv", row.names = 1)

# Remove all NAs
Kuster = Kuster[!is.na(Kuster$ProteinHalfLife),]
Mathieson = Mathieson[(rowSums(is.na(Mathieson[,2:5]))<4),]

# Merge Half-Life datasets
ENScons = unique(c(rownames(Mathieson),as.character(Kuster$EnsemblGeneID)))
halflife = data.frame(row.names = ENScons)
halflife$hela = sapply(ENScons,function(x) if(x %in% as.character(Kuster$EnsemblGeneID)){mean(Kuster[Kuster$EnsemblGeneID==x,"ProteinHalfLife"])}else{NA} )
halflife$Bcells = sapply(ENScons,function(x) if(x %in% rownames(Mathieson)){Mathieson[x,"BCells_ProHL"]}else{NA})
halflife$NKcells = sapply(ENScons,function(x) if(x %in% rownames(Mathieson)){Mathieson[x,"NKCells_ProHL"]}else{NA})
halflife$Hepatocytes = sapply(ENScons,function(x) if(x %in% rownames(Mathieson)){Mathieson[x,"Hepatocytes_ProHL"]}else{NA})
halflife$Monocytes = sapply(ENScons,function(x) if(x %in% rownames(Mathieson)){Mathieson[x,"Monocytes_ProHL"]}else{NA})
cor(halflife, method = "spearman", use = "pairwise.complete.obs")
write.csv(halflife,"data/ProteinHalfLife_merged.csv")

# For each tissue, check if there are differences in half-life
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t)&(TSprot$gene %in% ENScons),]
  dataset_temp = data.frame(row.names = seq(nrow(TStemp)*ncol(halflife)))
  dataset_temp$gene = rep(as.character(TStemp$gene),ncol(halflife))
  dataset_temp$halflife = as.numeric(unlist(halflife[as.character(TStemp$gene),]))
  dataset_temp$tissue = t
  dataset_temp$updown = rep(as.character(TStemp$updown),ncol(halflife))
  dataset_temp$HLcells = as.character(sapply(colnames(halflife), rep, nrow(TStemp)))
  dataset = rbind(dataset,dataset_temp)
}

# Plot
pdf("plots/ProtHalfLife_ALL.pdf",height=17,width =15)
diffexp = compare_means(halflife~updown, data=dataset, group.by = c("tissue","HLcells"), method="wilcox.test")
print(ggplot(dataset, aes(x=tissue, y=halflife, fill=updown)) +
        facet_wrap(. ~ HLcells, ncol = 1, scales = "free") +
        geom_boxplot(position=position_dodge(1)) +
        scale_fill_manual(values=c("blue", "red")) +
        labs(x="Tissue", y = "Protein Half-Life") +
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) +
        stat_compare_means(aes(label = ..p.adj..)) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
dev.off()

# Get average over all 5 HL datasets
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t)&(TSprot$gene %in% ENScons),]
  dataset_temp = data.frame(row.names = seq(nrow(TStemp)))
  dataset_temp$gene = as.character(TStemp$gene)
  dataset_temp$halflife = rowMeans(halflife[as.character(TStemp$gene),],na.rm=T)
  dataset_temp$tissue = t
  dataset_temp$updown = as.character(TStemp$updown)
  dataset = rbind(dataset,dataset_temp)
}

# Plot
pdf("plots/ProtHalfLife_avg_ALL.pdf",height=6,width =18)
diffexp = compare_means(halflife~updown, data=dataset, group.by = c("tissue"), method="wilcox.test")
print(ggplot(dataset, aes(x=tissue, y=halflife, fill=updown)) +
        geom_boxplot(position=position_dodge(1)) +
        labs(x="Tissue", y = "Protein Half-Life") +
        scale_fill_manual(values=c("blue", "red")) +
        stat_summary(position=position_dodge(1),fun.y=median, geom="point", shape=23, size=2) +
        stat_compare_means(aes(label = ..p.adj..)) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
dev.off()