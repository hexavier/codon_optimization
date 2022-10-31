library(ggplot2)
library(ggpubr)

# Load data
dataset=c()
files <- list.files(path="data/facs/primary", pattern="*.csv", full.names=TRUE, recursive=FALSE)
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
dataset$cell = sapply(dataset$cell, function(x) if(x=="l"){"Lung"}else 
      if(x=="r"){"Renal"})

# Compute ratio
dataset["ratio"] = log2(dataset$eGFP / dataset$mCherry)
dataset[is.infinite(dataset$ratio),"ratio"] <- NA

# Plot
ggplot(dataset, aes(x=plasmid,y=ratio,fill=cell)) +
  facet_grid(rep~., scales = "free_y") +
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
ggsave("plots/violin_facs_primary.pdf",width = 8, height = 4)

ggplot(dataset, aes(x=cell,y=ratio,fill=rep)) +
  facet_grid(~ plasmid) +
  geom_hline(yintercept=0) +
  geom_violin() +
  stat_summary(fun=median, geom="point", shape=23, size=2, position = position_dodge(width=0.9)) +
  stat_summary(aes(label=after_stat(y), y = stage(ratio, after_stat = 15)), 
               fun=length, geom="text", size=4, col = "black", vjust=-0.5,
               position = position_dodge(0.9)) +
  theme_classic() +
  labs(y = "log2(eGFP/mCherry)")
