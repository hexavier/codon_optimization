library(ggplot2)
library(ggpubr)
library(viridis)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## Load data
rnaseq = read.delim("data/RNAseq_Wang2020.txt", sep="\t", row.names = 1) # FPKM
rnaseq[rnaseq==0] <- NA
riboseq = read.delim("data/Ribo-seq_Wang2020.txt", sep="\t", row.names = 1) # FPKM
riboseq[riboseq==0] <- NA
ptr = read.delim("data/PTR_Eraslan2019.tsv", sep="\t", row.names = 2)[,-c(1,2,3)]

## Compute translational efficiency Wang 2020
tissues = sapply(colnames(rnaseq), function(x) strsplit(x,"_")[[1]][1])
medianRNA = sapply(unique(tissues),function(x) apply(log10(rnaseq[,tissues %in% x]),1,median, na.rm=T))
medianRibo = sapply(unique(tissues),function(x) apply(log10(riboseq[,tissues %in% x]),1,median, na.rm=T))
TE = medianRibo - medianRNA
# Clean NA and Inf
TE[is.na(TE)] <- NA
TE[is.infinite(TE)] <- NA
write.csv(TE,"data/TE_Wang2020.csv")

## Compare PTR vs TE
ENScons = rownames(TE)[rownames(TE) %in% rownames(ptr)]
tissues = colnames(ptr)[tolower(colnames(ptr)) %in% colnames(TE)]
dataset = c()
for (t in tissues){
  dataset_temp = data.frame(row.names = 1:length(ENScons))
  dataset_temp$TE = as.numeric(TE[ENScons,tolower(t)])
  dataset_temp$PTR = as.numeric(ptr[ENScons,t])
  dataset_temp$gene = ENScons
  dataset_temp$tissue = t
  # Merge
  dataset = rbind(dataset,dataset_temp)
}
dataset = dataset[rowSums(is.na(dataset[,c("TE","PTR")]))==0,]
dataset$density <- get_density(dataset$PTR, dataset$TE, n = 100)

# Plot
pdf("plots/TE-Wang2020_vs_PTR-Wang2019.pdf",height=7,width =3.8)
# Print correlations
print(ggplot(dataset, aes(x=PTR, y=TE, color = density)) +
        facet_wrap(. ~ tissue, ncol = 1, scales = "free") +
        geom_point(size=1, alpha=1,shape=19, stroke=0,na.rm=F) + 
        scale_color_viridis() +
        stat_cor(method = "spearman") +
        theme_classic())
dev.off()
