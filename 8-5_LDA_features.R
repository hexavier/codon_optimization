library(ggplot2)
library(ggpubr)
library(gplots)

# Load data
lda = t(read.csv("results/LDA_classes.csv", row.names = 1))
codusnormratios = read.csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
codus_UPDOWN = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)[,3:62]

# Define structure
tissues = rownames(lda)[rownames(lda) %in% rownames(codus_UPDOWN)]
codons = colnames(codus_UPDOWN); codons = codons[-4]
dataset = c()
for (t in tissues){
  dataset_temp = data.frame(row.names = 1:length(codons))
  dataset_temp$lda = as.numeric(lda[t,codons])
  dataset_temp$codus_UPDOWN = as.numeric(codus_UPDOWN[t,codons])
  dataset_temp$codusnormratios = as.numeric(codusnormratios[t,codons])
  dataset_temp$tissue = t
  dataset_temp$codon = codons
  # Merge
  dataset = rbind(dataset,dataset_temp)
}

# Plot LDA classes
heatmap.2(t(lda),symm=F, margins = c(6,5), col = bluered, na.color = "grey")


pdf("plots/RiboSeqLDA_vs_TisOpt.pdf",height=8,width =12)

# Print correlations
print(ggplot(dataset, aes(x=lda, y=codusnormratios)) +
        facet_wrap(. ~ tissue, ncol = 3, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=abs(lda), y=codus_UPDOWN)) +
        facet_wrap(. ~ tissue, ncol = 3, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
dev.off()

# Correlation table
cor(t(rbind(lda[,codons],codusnormratios[,codons])),method="spearman",use = "pairwise")[(nrow(lda)+1):(nrow(lda)+nrow(codusnormratios)),1:nrow(lda)]
