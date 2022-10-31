library(ggplot2)
library(ggpubr)
library(gplots)

# Load data
all_ratios = read.csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
all_aucs = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)
mrna_ratios = read.csv("results/all-mRNA_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
mrna_aucs = read.csv("results/mean_mRNA-all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)
prot_ratios = read.csv("results/all-prot_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
prot_aucs = read.csv("results/mean_prot-all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)
fcglob0_ratios = read.csv("results/all-FCglob0_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
fcglob0_aucs = read.csv("results/mean_FCglob0-all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)

# ALL dataset
outratios <- heatmap.2(t(all_ratios),col=bluered, symm=F, margins = c(6,5), na.color = "grey", trace="none")
tissueorder = rownames(all_ratios)[outratios$colInd]
codonorder = rev(colnames(all_ratios)[outratios$rowInd])

# Plot mRNA with same order as PTR
mrna_ordered = mrna_ratios[tissueorder[tissueorder %in% rownames(mrna_ratios)],codonorder[codonorder %in% colnames(mrna_ratios)]]
outratios <- heatmap.2(t(mrna_ordered), col=bluered, Rowv = F, Colv = F, symm=F, margins = c(6,5), na.color = "grey", trace="none")
aucs = mrna_aucs[tissueorder[tissueorder %in% rownames(mrna_ratios)],1]
barplot(aucs,ylim=c(0,1)); abline(h=0.5)

# Plot prot with same order as PTR
prot_ordered = prot_ratios[tissueorder[tissueorder %in% rownames(prot_ratios)],codonorder[codonorder %in% colnames(prot_ratios)]]
outratios <- heatmap.2(t(prot_ordered), col=bluered, Rowv = F, Colv = F, symm=F, margins = c(6,5), na.color = "grey", trace="none")
aucs = prot_aucs[tissueorder[tissueorder %in% rownames(prot_ratios)],1]
barplot(aucs,ylim=c(0,1)); abline(h=0.5)

# Plot FCglob0 with same order as PTR
prot_ordered = fcglob0_ratios[tissueorder[tissueorder %in% rownames(fcglob0_ratios)],codonorder[codonorder %in% colnames(fcglob0_ratios)]]
outratios <- heatmap.2(t(prot_ordered), col=bluered, Rowv = F, Colv = F, symm=F, margins = c(6,5), na.color = "grey", trace="none")
aucs = fcglob0_aucs[tissueorder[tissueorder %in% rownames(fcglob0_ratios)],1]
barplot(aucs,ylim=c(0,1)); abline(h=0.5)