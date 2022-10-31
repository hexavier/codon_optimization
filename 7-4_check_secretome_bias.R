library(ggplot2)
library(ggpubr)

# Load data
TSprot = read.csv("results/TSall_PTR-UP-FCavg2-FCglob1.tsv", sep="\t")
# 1708 genes secretome: https://www.proteinatlas.org/humanproteome/tissue/secretome
HPA_predictions = read.csv("data/HPA_secretome.tsv", sep="\t")

# Get list of secreted proteins
secretome = unique(as.character(HPA_predictions$Ensembl))

# For each tissue, check secretion
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t),]
  dataset_temp = data.frame(row.names = seq(4))
  dataset_temp$updown = c("UP","UP","DOWN","DOWN")
  dataset_temp$secreted = c(T,F,T,F)
  dataset_temp$tissue = t
  dataset_temp$number = apply(dataset_temp,1,function(x) sum((TStemp[TStemp$updown==x["updown"],"gene"] %in% secretome)==x["secreted"]))
  dataset = rbind(dataset,dataset_temp)
}

# Plot
ggplot(dataset, aes(x=tissue, y=number, fill=updown)) +
  facet_wrap(. ~ secreted, ncol = 1, scales = "free") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("blue", "red")) +
  labs(x="Tissue", y = "# Proteins") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

### SECRETED TO BLOOD ###
# 729 genes secreted to blood: https://www.proteinatlas.org/humanproteome/blood/secretome
HPA_toBlood = read.csv("data/HPA_toBlood.tsv", sep="\t")
toBlood = unique(as.character(HPA_toBlood$Ensembl))
# Build dataset
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t),]
  dataset_temp = data.frame(row.names = seq(4))
  dataset_temp$updown = c("UP","UP","DOWN","DOWN")
  dataset_temp$secreted = c(T,F,T,F)
  dataset_temp$tissue = t
  dataset_temp$number = apply(dataset_temp,1,function(x) sum((TStemp[TStemp$updown==x["updown"],"gene"] %in% toBlood)==x["secreted"]))
  dataset = rbind(dataset,dataset_temp)
}
# Plot
ggplot(dataset, aes(x=tissue, y=number, fill=updown)) +
  facet_wrap(. ~ secreted, ncol = 1, scales = "free") +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("blue", "red")) +
  labs(x="Tissue", y = "# Proteins") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))