library("FactoMineR")
library("factoextra")
library(sva)

# Load data
codocc = t(read.csv("data/codon_occupancies.csv", row.names = 1))
metadata = read.delim("data/riboseq_datasets.tsv", sep="\t", row.names = 1)

# ComBat
map_batch = c(1,2,3,4,5); names(map_batch) = c("Zou 2019","Gonzalez 2014","Loayza-Puch 2016","van Heesch 2019","Wang 2020")
batch = sapply(colnames(codocc),function(x) map_batch[metadata[x,"study"]])
combat_codocc = ComBat(dat=codocc, batch=batch, prior.plots =TRUE)

# Remove NAs and stop codons
dataset = combat_codocc
dataset[dataset==0] <- NA
dataset = dataset[-grep("TAA|TAG|TGA",rownames(dataset)),]

# Correspondence analysis
res.ca <- CA(dataset, graph = FALSE)
# Symetric plot
fviz_ca_biplot(res.ca, repel = TRUE)


### Median per tissue
tissues = unique(metadata$tissue)
medians = data.frame(row.names = rownames(dataset))
for (t in tissues){
  samples = colnames(dataset)[colnames(dataset) %in% rownames(metadata)[metadata$tissue %in% t]]
  medians[,t] = apply(dataset[,samples],1,median,na.rm=T)
}

# Correspondence analysis
res.ca <- CA(medians, graph = FALSE)
# Symetric plot
fviz_ca_biplot(res.ca, repel = TRUE, col.col = "black", col.row= "grey")
