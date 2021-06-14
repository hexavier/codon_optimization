library(ggplot2)
library(ggpubr)
library(gplots)
library(tAI)

extract_cod <- function (trnas, anticod){
  output = data.frame(row.names = anticod)
  trnas_acod = sapply(rownames(trnas), function(x) substr(x,nchar(x)-2,nchar(x)))
  for (s in colnames(trnas)){
    output[,s] = sapply(anticod, function(x) if(any(trnas_acod==x)){mean(trnas[trnas_acod==x,s])}else{0})
  }
  return(output)
}

# Load data
codocc = t(read.csv("data/codon_occupancies.csv", row.names = 1))
metadata = read.delim("data/riboseq_datasets.tsv", sep="\t", row.names = 1)
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)

# Remove NAs and stop codons
codocc[codocc==0] <- NA
codocc = codocc[-grep("TAA|TAG|TGA",rownames(codocc)),]
codocc = log2(codocc)

# Plot raw
heatmap.2(codocc,symm=F, margins = c(6,5), col = bluered, na.color = "grey")

# Check correlation with tAI
copy_numbers = read.delim("data/tRNA_copy_numbers.tsv", sep="\t", row.names = 1)
anticodon = extract_cod(copy_numbers,codons$ANTICODON)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
tAI = get.ws(tRNA=anticodon[,1], s=initial_s, sking=0)
codonocc_median = apply(extract_cod(codocc,rownames(codons)[!(codons$AA %in% c("Stop","Met"))])[,c(1,2)],1,median,na.rm=T)
cor(codonocc_median,tAI, method="spearman")

# Median per tissue
tissues = unique(metadata$tissue)
medians = data.frame(row.names = rownames(codocc))
for (t in tissues){
  samples = colnames(codocc)[colnames(codocc) %in% rownames(metadata)[metadata$tissue %in% t]]
  medians[,t] = apply(codocc[,samples],1,median,na.rm=T)
}

# Plot raw
heatmap.2(as.matrix(medians),symm=F, margins = c(6,5), col = bluered, na.color = "grey")

## Check correlation medians with tAI
codons = read.csv("data/codons_table.tab", sep="\t", row.names = 1)
copy_numbers = read.delim("data/tRNA_copy_numbers.tsv", sep="\t", row.names = 1)
anticodon = extract_cod(copy_numbers,codons$ANTICODON)
initial_s = c(0, 0, 0, 0, 0.5, 0.5, 0.75, 0.5, 0.5)
tAI = get.ws(tRNA=anticodon[,1], s=initial_s, sking=0)
# Ordered occ
codonord = extract_cod(medians,rownames(codons)[!(codons$AA %in% c("Stop","Met"))])
# Correlate
cor(cbind(codonord,tAI), method = "spearman", use = "pairwise")[1:ncol(codonord),(ncol(codonord)+1)]

## Check correlation with tissue codon usage
codusWang = read.csv("results/PTR_tissue_codus.csv", row.names = 1)
codusord = extract_cod(codusWang,rownames(codons)[!(codons$AA %in% c("Stop","Met"))])
cor(cbind(codonord,codusord), method = "spearman", use = "pairwise")[(ncol(codonord)+1):(ncol(codonord)+ncol(codusord)),1:ncol(codonord)]


### COMPARE TisOpt ###
# Load data
codusnormratios = read.csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", row.names = 1)
codus_UPDOWN = read.csv("results/mean_all_randomforest_codusnorm_to_UPDOWNproteins.csv", row.names = 1)[,3:62]
medians = t(medians)

# Create structure
tissues = tissues[tissues %in% rownames(codus_UPDOWN)]
codons = colnames(codus_UPDOWN); codons = codons[-4]
dataset = c()
for (t in tissues){
  dataset_temp = data.frame(row.names = 1:length(codons))
  dataset_temp$codonocc = as.numeric(medians[t,codons])
  dataset_temp$codus_UPDOWN = as.numeric(codus_UPDOWN[t,codons])
  dataset_temp$codusnormratios = as.numeric(codusnormratios[t,codons])
  dataset_temp$tissue = t
  dataset_temp$codon = codons
  # Merge
  dataset = rbind(dataset,dataset_temp)
}

# Plot
pdf("plots/codonocc_vs_TisOpt.pdf",height=8,width =12)

# Print correlations
print(ggplot(dataset, aes(x=codonocc, y=codusnormratios)) +
        facet_wrap(. ~ tissue, ncol = 3, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
print(ggplot(dataset, aes(x=abs(codonocc), y=codus_UPDOWN)) +
        facet_wrap(. ~ tissue, ncol = 3, scales = "free") +
        geom_point(size=2,alpha=0.7,na.rm=F) + 
        stat_cor(method = "spearman") +
        theme_classic())
dev.off()
