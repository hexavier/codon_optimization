library(ggplot2)
library(ggpubr)
library(gplots)

# Load data
all_ratios = read.csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
all_aucs = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)
gtex_ratios = read.csv("results/ptr-gtex_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
gtex_aucs = read.csv("results/mean_ptr-gtex_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)
eraslan_ratios = read.csv("results/ptr-eraslan_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
eraslan_aucs = read.csv("results/mean_ptr-eraslan_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)

# ALL dataset
outratios <- heatmap.2(t(all_ratios),col=bluered, symm=F, margins = c(6,5), na.color = "grey", trace="none")
tissueorder = rownames(all_ratios)[outratios$colInd]
codonorder = rev(colnames(all_ratios)[outratios$rowInd])

# Plot GTEx with same order as ALL
gtex_ordered = gtex_ratios[tissueorder[tissueorder %in% rownames(gtex_ratios)],codonorder[codonorder %in% colnames(gtex_ratios)]]
outratios <- heatmap.2(t(gtex_ordered), col=bluered, Rowv = F, Colv = F, symm=F, margins = c(6,5), na.color = "grey", trace="none")
aucs = gtex_aucs[tissueorder[tissueorder %in% rownames(gtex_ratios)],1]
barplot(aucs,ylim=c(0,1)); abline(h=0.5)

# Plot Eraslan with same order as ALL
eraslan_ordered = eraslan_ratios[tissueorder[tissueorder %in% rownames(eraslan_ratios)],codonorder[codonorder %in% colnames(eraslan_ratios)]]
outratios <- heatmap.2(t(eraslan_ordered), col=bluered, Rowv = F, Colv = F, symm=F, margins = c(6,5), na.color = "grey", trace="none")
aucs = eraslan_aucs[tissueorder[tissueorder %in% rownames(eraslan_ratios)],1]
barplot(aucs,ylim=c(0,1)); abline(h=0.5)
