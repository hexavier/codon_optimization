library(gplots)

# Load data
coduspairratios = read.csv("results/all_CodPairs_ratios_UPDOWNproteins.csv", row.names = 1)

# Remove infinites and NA
coduspairratios[is.infinite(as.matrix(coduspairratios))] <- NA
coduspairratios = coduspairratios[colSums(is.na(coduspairratios))<nrow(coduspairratios)]
coduspairratios[is.na(as.matrix(coduspairratios))] <- 0

# Keep only top 100 extreme bicodons
ranks = order(rowSums(apply(abs(as.matrix(coduspairratios)),1,rank)),decreasing = T)[1:100]
heatmap.2(as.matrix(coduspairratios)[,ranks], col="bluered", margins = c(6,5), na.color = "grey")

# Get correlation table
cortable = cor(t(coduspairratios), method = "spearman", use = "pairwise")
heatmap.2(cortable, col="bluered", margins = c(7,7))

#### Compare with individual codons ####
library(ggplot2)
library(ggpubr)
# Load data
codusnormratios = read.csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
coduspairratios = read.csv("results/all_CodPairs_ratios_UPDOWNproteins.csv", row.names = 1)
# Remove infinites and NA
coduspairratios[is.infinite(as.matrix(coduspairratios))] <- NA
coduspairratios = coduspairratios[colSums(is.na(coduspairratios))<nrow(coduspairratios)]

# Compute expected pairs based on addition of individual codons
expectedpairs = data.frame(row.names = rownames(coduspairratios))
for (p in colnames(coduspairratios)){
  codons = substring(toupper(p), c(1,4), c(3,6))
  if (all(codons %in% colnames(codusnormratios))){
    expectedpairs[,p] = as.numeric(rowSums(codusnormratios[rownames(coduspairratios),codons]))
  }
}

# Correlate in each tissue
dataset = c()
for (t in rownames(expectedpairs)){
  dataset_temp = data.frame(row.names = 1:ncol(expectedpairs))
  dataset_temp$expected = as.numeric(expectedpairs[t,])
  dataset_temp$observed = as.numeric(coduspairratios[t,colnames(expectedpairs)])
  dataset_temp$codpair = as.character(colnames(expectedpairs))
  dataset_temp$tissue = t
  dataset = rbind(dataset,dataset_temp)
}

### Identify points deviating from linear model ###
codoncooksd = data.frame(row.names = colnames(expectedpairs))
codonresid = data.frame(row.names = colnames(expectedpairs))
lmsummary = data.frame(row.names = c("R2","RSE","pvalue"))
for (t in rownames(expectedpairs)){
  lm.r = lm(observed ~ expected, data = dataset[dataset$tissue==t,])
  codoncooksd[as.character(dataset[names(cooks.distance(lm.r)),"codpair"]),t] = as.numeric(cooks.distance(lm.r))
  codonresid[as.character(dataset[names(rstandard(lm.r)),"codpair"]),t] = as.numeric(rstandard(lm.r))
  lmsummary["RSE",t] = summary(lm.r)[[6]]
  lmsummary["R2",t] = summary(lm.r)[[9]]
  fstat <- summary(lm.r)$fstatistic
  lmsummary["pvalue",t] = as.numeric(pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE))
  dataset[names(rstandard(lm.r)),"bias"] = as.numeric(rstandard(lm.r))
}
write.csv(codoncooksd, "results/TisOpt_CodPairs_tissue_CooksD.csv", quote=F)
write.csv(codonresid, "results/TisOpt_CodPairs_tissue_rstandard.csv", quote=F)
write.csv(lmsummary, "results/TisOpt_CodPairs_tissue_regression.csv", quote=F)


# Show top outliers

##Plot outliers
ranks = order(rowSums(apply(t(abs(codonresid)),1,rank)),decreasing = T)[1:100]
heatmap.2(t(codonresid)[,ranks], col="bluered", margins = c(4,8), na.color = "grey")
heatmap.2(cor(codonresid, method = "spearman", use="pairwise"), col="bluered", margins = c(8,8), na.color = "grey")

pdf("plots/correlation_CodPairs_ObsVSExp.pdf",height=15,width =20)
print(ggplot(dataset, aes(x=expected, y=observed)) +
        facet_wrap(. ~ tissue, ncol = 6, scales = "free") +
        geom_point(size=2,alpha=0.1,na.rm=F) + 
        geom_smooth(method=lm) +
        geom_text(aes(label=ifelse(abs(bias)>4,as.character(codpair),'')),hjust=0,vjust=0,size=3) +
        stat_cor(method = "spearman") +
        theme_classic())
dev.off()

### Check codpair counts of outliers
upcodpaircounts = read.csv("results/all_CodPairs_counts_UPproteins.csv", row.names = 1)
downcodpaircounts = read.csv("results/all_CodPairs_counts_DOWNproteins.csv", row.names = 1)
# Add counts to dataset
updataset = dataset[,c("bias","codpair","tissue")]; downdataset = dataset[,c("bias","codpair","tissue")]
updataset["updown"] = "UP"; downdataset["updown"] = "DOWN";
for (t in rownames(upcodpaircounts)){
  tissueidx = which(updataset$tissue==t)
  updataset[tissueidx,"counts"] = as.numeric(upcodpaircounts[t,as.character(updataset[tissueidx,"codpair"])])
  downdataset[tissueidx,"counts"] = as.numeric(downcodpaircounts[t,as.character(downdataset[tissueidx,"codpair"])])
}
updataset["r_bins"] = cut(updataset$bias,breaks=9); downdataset["r_bins"] = cut(downdataset$bias,breaks=9)
alldataset = rbind(updataset,downdataset)

# Plot
ggplot(alldataset[!is.na(alldataset$r_bins),], aes(x=r_bins, y=counts, fill=updown)) +
  geom_boxplot(position=position_dodge(1)) + 
  scale_fill_manual(values=c("blue","red")) +
  scale_y_continuous(trans='log10') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))
cor(abs(updataset$bias),updataset$counts,method="spearman", use="pairwise")
cor(abs(downdataset$bias),downdataset$counts,method="spearman", use="pairwise")
