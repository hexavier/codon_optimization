library(ggplot2)
library(ggpubr)
library(reshape)
library(ggrepel)
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

# Compute ratio
dataset["ratio"] = log2(dataset$eGFP / dataset$mCherry)

# Write full names
dataset$plasmid = sapply(dataset$plasmid, function(x) if(x=="74"){"mCherry_L-eGFP_K-1"}else 
  if(x=="75"){"mCherry_K-eGFP_L-1"}else 
    if(x=="76"){"mCherry_L-eGFP_K-2"} else 
      if(x=="77"){"mCherry_K-eGFP_L-2"})
dataset$plasmid = factor(as.character(dataset$plasmid), levels=c("mCherry_L-eGFP_K-1","mCherry_L-eGFP_K-2","mCherry_K-eGFP_L-1","mCherry_K-eGFP_L-2"))
dataset$cell = sapply(dataset$cell, function(x) if(x=="h"){"HEK293T"}else 
  if(x=="a"){"A549"}else 
    if(x=="l"){"Lung"}else 
      if(x=="r"){"Renal"}else 
        if(x=="z"){"Optimem"})

# Compute mean and SD
mediandf1 = data.frame()
mediandf2 = data.frame()
n=1
for (ce in unique(dataset$cell)){
  for (r in unique(dataset[dataset$cell==ce,"rep"])){
    tempdf = dataset[(dataset$cell==ce)&(dataset$rep==r),]
    mediandf1[n,c("cell","rep","xmedian","xup","xdown",
                  "ymedian","yup","ydown")] = c(ce,sprintf("%s-%s",ce,r),
                                                median(tempdf[tempdf$plasmid=="mCherry_K-eGFP_L-1","ratio"],na.rm=T),
                                                quantile(tempdf[tempdf$plasmid=="mCherry_K-eGFP_L-1","ratio"],probs=0.6,na.rm=T)[1],
                                                quantile(tempdf[tempdf$plasmid=="mCherry_K-eGFP_L-1","ratio"],probs=0.4,na.rm=T)[1],
                                                median(tempdf[tempdf$plasmid=="mCherry_L-eGFP_K-1","ratio"],na.rm=T),
                                                quantile(tempdf[tempdf$plasmid=="mCherry_L-eGFP_K-1","ratio"],probs=0.6,na.rm=T)[1],
                                                quantile(tempdf[tempdf$plasmid=="mCherry_L-eGFP_K-1","ratio"],probs=0.4,na.rm=T)[1])
    mediandf2[n,c("cell","rep","xmedian","xup","xdown",
                  "ymedian","yup","ydown")] = c(ce,sprintf("%s-%s",ce,r),
                                                median(tempdf[tempdf$plasmid=="mCherry_K-eGFP_L-2","ratio"],na.rm=T),
                                                quantile(tempdf[tempdf$plasmid=="mCherry_K-eGFP_L-2","ratio"],probs=0.6,na.rm=T)[1],
                                                quantile(tempdf[tempdf$plasmid=="mCherry_K-eGFP_L-2","ratio"],probs=0.4,na.rm=T)[1],
                                                median(tempdf[tempdf$plasmid=="mCherry_L-eGFP_K-2","ratio"],na.rm=T),
                                                quantile(tempdf[tempdf$plasmid=="mCherry_L-eGFP_K-2","ratio"],probs=0.6,na.rm=T)[1],
                                                quantile(tempdf[tempdf$plasmid=="mCherry_L-eGFP_K-2","ratio"],probs=0.4,na.rm=T)[1])
    n=n+1
  }
}
mediandf1[,c("xmedian","xup","xdown","ymedian","yup","ydown")] = apply(mediandf1[,c("xmedian","xup","xdown","ymedian","yup","ydown")],2,as.numeric)
mediandf2[,c("xmedian","xup","xdown","ymedian","yup","ydown")] = apply(mediandf2[,c("xmedian","xup","xdown","ymedian","yup","ydown")],2,as.numeric)

# Plot
p1 <-ggplot(mediandf1, aes(x=xmedian,y=ymedian,label=rep,color=cell)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ydown, ymax = yup, width=.1)) + 
  geom_errorbarh(aes(xmin = xdown, xmax = xup, height=.1)) +
  scale_color_manual(values=c("#304d63","#8fb9aa")) +
  geom_text_repel() +
  theme_classic() +
  labs(title="Plasmids Set 1", subtitle="log2(eGFP/mCherry)",x="mCherry_K-eGFP_L", y="mCherry_L-eGFP_K")

p2 <- ggplot(mediandf2, aes(x=xmedian,y=ymedian,label=rep,color=cell)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ydown, ymax = yup, width=.1)) + 
  geom_errorbarh(aes(xmin = xdown, xmax = xup, height=.1)) +
  scale_color_manual(values=c("#304d63","#8fb9aa")) +
  geom_text_repel() +
  theme_classic() +
  labs(title="Plasmids Set 2", subtitle="log2(eGFP/mCherry)",x="mCherry_K-eGFP_L", y="mCherry_L-eGFP_K")
p1 | p2
ggsave("plots/scatter_facs_celllines.pdf",width = 10, height = 4)
