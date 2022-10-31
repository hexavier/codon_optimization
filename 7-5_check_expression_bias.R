library(ggplot2)
library(ggpubr)

# Load data
TSprot = read.csv("results/TSall_PTR-UP-FCavg2-FCglob1.tsv", sep="\t")

### mRNA ###
# Load Eraslan
mrnaeraslan = read.csv("data/mRNAseq_FPKM_Eraslan2019.csv",row.names = 2)[,5:33]

# Load GTEx
mrnagtex = read.csv("data/mRNAseq_GTEx_Jiang2020.csv",row.names = 1,skip = 1)
mapping = read.csv("data/GTEx_samples.csv")
# Remove TPM<10 (same as Eraslan)
mrnagtex[mrnagtex<log2(10)] = NA; mrnagtex = mrnagtex[rowSums(!is.na(mrnagtex))>0,]
mrna_samples = sapply(colnames(mrnagtex),function(x) paste0(unlist(strsplit(x,"\\."))[1:3],collapse="."))
colnames(mrnagtex) = mrna_samples
tissues = as.character(sapply(colnames(mrnagtex), function(x) mapping[mapping$Sample.ID==x,"Tissue"]))
mrna_med = sapply(unique(tissues), function(x) apply(mrnagtex[,tissues %in% x, drop=F],1,median, na.rm=T))

# For each tissue, check if there are differences in expression
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t),]
  dataset_temp = data.frame(row.names = seq(nrow(TStemp)))
  dataset_temp$gene = as.character(TStemp$gene)
  if (t %in% colnames(mrnaeraslan)){
    dataset_temp$Eraslan = as.numeric(unlist(sapply(as.character(TStemp$gene),function(x)if(x %in% rownames(mrnaeraslan)){mrnaeraslan[x,t]}else{NA})))
  }else{
    dataset_temp$Eraslan = NA
  }
  if (t %in% colnames(mrna_med)){
    dataset_temp$GTEx = as.numeric(unlist(sapply(as.character(TStemp$gene),function(x)if(x %in% rownames(mrna_med)){mrna_med[x,t]}else{NA})))
  }else{
    dataset_temp$GTEx = NA
  }
  dataset_temp$tissue = t
  dataset_temp$updown = as.character(TStemp$updown)
  dataset = rbind(dataset,dataset_temp)
}

# Plot
pdf("plots/expression_updown_proteins.pdf",height=6,width =16)
print(ggplot(dataset, aes(x=tissue, y=GTEx, fill=updown)) +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values=c("blue", "red")) +
  labs(x="Tissue", y = "mRNA expresion") +
  stat_summary(position=position_dodge(1),fun=median, geom="point", shape=23, size=2) +
  stat_compare_means(aes(label = ..p.adj..)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))
print(ggplot(dataset, aes(x=tissue, y=Eraslan, fill=updown)) +
  geom_boxplot(position=position_dodge(1)) +
  scale_fill_manual(values=c("blue", "red")) +
  labs(x="Tissue", y = "mRNA expresion") +
  stat_summary(position=position_dodge(1),fun=median, geom="point", shape=23, size=2) +
  stat_compare_means(aes(label = ..p.adj..)) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))

### PROTEIN ###
# Load Eraslan (GTEx is relative quantification, not useful here)
proteraslan = read.csv("data/proteome_iBAQ_Eraslan2019.csv", row.names = 2)[,5:33]

# For each tissue, check if there are differences in expression
dataset = c()
for(t in unique(TSprot$Tissue)){
  TStemp = TSprot[(TSprot$Tissue==t),]
  dataset_temp = data.frame(row.names = seq(nrow(TStemp)))
  dataset_temp$gene = as.character(TStemp$gene)
  if (t %in% colnames(proteraslan)){
    dataset_temp$Eraslan = as.numeric(unlist(sapply(as.character(TStemp$gene),function(x)if(x %in% rownames(proteraslan)){proteraslan[x,t]}else{NA})))
  }else{
    dataset_temp$Eraslan = NA
  }
  dataset_temp$tissue = t
  dataset_temp$updown = as.character(TStemp$updown)
  dataset = rbind(dataset,dataset_temp)
}

# Plot
print(ggplot(dataset, aes(x=tissue, y=Eraslan, fill=updown)) +
        geom_boxplot(position=position_dodge(1)) +
        scale_fill_manual(values=c("blue", "red")) +
        labs(x="Tissue", y = "Protein expresion") +
        stat_summary(position=position_dodge(1),fun=median, geom="point", shape=23, size=2) +
        stat_compare_means(aes(label = ..p.adj..)) +
        theme_classic() + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)))

dev.off()

