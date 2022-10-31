library(ggplot2)
library(ggpubr)
library(patchwork)

# Load data
dataset=c()
files <- list.files(path="data/facs/celllines", pattern="*.csv", full.names=TRUE, recursive=FALSE)
for (f in files){
  rawfile = read.delim(f,header=F,skip=1,sep="")
  tempfile = as.data.frame(t(apply(rawfile,1,function(x) c(as.numeric(paste0(strsplit(x,",")[[1]][1:2],collapse=".")),
                                                           as.numeric(paste0(strsplit(x,",")[[1]][3:4],collapse="."))))))
  colnames(tempfile) = c("eGFP","mCherry")
  tempfile["file"] = f
  tempfile["cell"] = substr(strsplit(f,"/|_")[[1]][6],1,1)
  tempfile["plasmid"] = substr(strsplit(f,"/|_")[[1]][6],2,3)
  tempfile["rep"] = substr(strsplit(f,"/|_")[[1]][7],1,4)
  dataset=rbind(dataset,tempfile)
}
# Write full names
dataset$plasmid = sapply(dataset$plasmid, function(x) if(x=="74"){"mCherry_L-eGFP_K-1"}else 
  if(x=="75"){"mCherry_K-eGFP_L-1"}else 
    if(x=="76"){"mCherry_L-eGFP_K-2"} else 
      if(x=="77"){"mCherry_K-eGFP_L-2"})
dataset$plasmid = factor(as.character(dataset$plasmid), levels=c("mCherry_L-eGFP_K-1","mCherry_L-eGFP_K-2","mCherry_K-eGFP_L-1","mCherry_K-eGFP_L-2"))
dataset$cell = sapply(dataset$cell, function(x) if(x=="h"){"HEK293T"}else 
  if(x=="a"){"A549"})

# Compute ratio
dataset["ratio"] = log2(dataset$eGFP / dataset$mCherry)
dataset = dataset[!is.na(dataset$ratio)&!is.infinite(dataset$ratio),]

# Plot
ggplot(dataset, aes(x=plasmid,y=ratio,fill=cell)) +
  facet_grid(rep~.)+
  geom_hline(yintercept=0) +
  geom_violin() +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  stat_summary(fun=median, geom="point", shape=23, size=2, position = position_dodge(width=0.9)) +
  stat_summary(aes(label=after_stat(y), y = stage(ratio, after_stat = 15)), 
               fun=length, geom="text", size=4, col = "black", vjust=-0.5,
               position = position_dodge(0.9)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_classic() +
  labs(y = "log2(eGFP/mCherry)")
ggsave("plots/violin_facs_cellines.pdf",width = 8, height = 5)

### MERGED PLOT ###
# To make sure that all replicates contribute equivalently, downsampling 1000 events per group
downdataset = c()
labels = sprintf("%s-%s-%s",dataset$cell,dataset$plasmid, dataset$rep)
for (g in unique(labels)){
  tempdf = dataset[labels==g,]
  downdataset = rbind(downdataset,tempdf[sample(rownames(tempdf),size=1000,replace=F),])
}
# Plot
p1<-ggplot(downdataset, aes(x=plasmid,y=ratio,fill=cell)) +
  geom_hline(yintercept=0) +
  geom_violin() +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  stat_summary(fun=median, geom="point", shape=23, size=2, position = position_dodge(width=0.9)) +
  stat_summary(aes(label=after_stat(y), y = stage(ratio, after_stat = 15)), 
               fun=length, geom="text", size=4, col = "black", vjust=-0.5,
               position = position_dodge(0.9)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  theme_classic() +
  theme(legend.position="top") +
  ylim(-10,15)+
  labs(y = "log2(eGFP/mCherry)")
p2<-ggplot(downdataset, aes(y=ratio,fill=cell)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  theme_classic() +
  ylim(-10,15)+
  theme(title = element_blank(),text=element_blank(),axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),legend.position="none")
p1 + p2 + plot_layout(widths = c(8, 1))
ggsave("plots/violin_facs_cellines_merged.pdf",width = 7, height = 3.5)
