library(ggplot2)
library(ggpubr)
library(gplots)

# Load data
codus_UPDOWN = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)[,3:62]
rcu = read.csv("results/mean_all_randomforest_RCU_to_UPDOWNproteins.csv", row.names = 1)[,3:62]
featTCGAsamples = read.csv("results/mean_randomforest_SDAs.csv", row.names = 1)[,3:62]
colnames(featTCGAsamples) = substr(colnames(featTCGAsamples),4,6)
ratiosSDA = read.csv("results/mean_RATIO_SDAs.csv", row.names = 1)
colnames(ratiosSDA) = substr(colnames(ratiosSDA),4,6)
codusnormratios = read.csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
RCUratios = read.csv("results/all_RCU_ratios_UPDOWNproteins.csv", row.names = 1)

# Mapping tissues
mapping = c(Tonsil="oral",Uterus="uterus",Prostate="prostate",
            Liver="liver",Colon="intestine",Brain="brain",
            Breast="breast",Esophagus="esophagus",Thyroid="thyroid",
            Stomach="stomach",Lung="lung",Skin="skin", Kidney="kidney",
            Urinarybladder="bladder")

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
  if(t %in% names(mapping)){
    dataset_temp$featTCGAsamples = as.numeric(featTCGAsamples[mapping[t],codons])
    dataset_temp$ratiosSDA = as.numeric(ratiosSDA[mapping[t],codons])
  }else{
    dataset_temp$featTCGAsamples = NA
    dataset_temp$ratiosSDA = NA
  }
  # Merge
  dataset = rbind(dataset,dataset_temp)
}

pdf("plots/ALL_codonSDA_vs_othermetrics.pdf",height=15,width =20)

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
print(ggplot(dataset, aes(x=codus_UPDOWN, y=featTCGAsamples)) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=featTCGAsamples, y=abs(ratiosSDA))) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=codus_UPDOWN, y=abs(ratiosSDA))) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=codusnormratios, y=ratiosSDA)) +
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
print(ggplot(dataset, aes(x=rcu, y=featTCGAsamples)) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=RCUratios, y=ratiosSDA)) +
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