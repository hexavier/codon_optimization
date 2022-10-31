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
mediandf = data.frame()
n=1
for (pl in unique(dataset$plasmid)){
  for (r in c("rep1","rep2","rep3")){
    tempdf = dataset[(dataset$plasmid==pl)&(dataset$rep==r),]
    mediandf[n,c("plasmid","rep","xmedian","xup","xdown",
                  "ymedian","yup","ydown")] = c(pl,sprintf("%s-%s",pl,r),
                                                median(tempdf[tempdf$cell=="HEK293T","ratio"],na.rm=T),
                                                quantile(tempdf[tempdf$cell=="HEK293T","ratio"],probs=0.6,na.rm=T)[1],
                                                quantile(tempdf[tempdf$cell=="HEK293T","ratio"],probs=0.4,na.rm=T)[1],
                                                median(tempdf[tempdf$cell=="A549","ratio"],na.rm=T),
                                                quantile(tempdf[tempdf$cell=="A549","ratio"],probs=0.6,na.rm=T)[1],
                                                quantile(tempdf[tempdf$cell=="A549","ratio"],probs=0.4,na.rm=T)[1])
    n=n+1
  }
}
mediandf[,c("xmedian","xup","xdown","ymedian","yup","ydown")] = apply(mediandf[,c("xmedian","xup","xdown","ymedian","yup","ydown")],2,as.numeric)

# Plot
ggplot(mediandf, aes(x=xmedian,y=ymedian,label=rep,color=plasmid)) +
  geom_point() + 
  geom_errorbar(aes(ymin = ydown, ymax = yup, width=.1)) + 
  geom_errorbarh(aes(xmin = xdown, xmax = xup, height=.1)) +
  scale_color_manual(values=c("#304d63", "#304d63", "#8fb9aa", "#8fb9aa")) +
  geom_text_repel(size=3) +
  theme_classic() +
  ylim(-10,10) +
  xlim(-10,10) +
  labs(title="log2(eGFP/mCherry)",x="HEK293T", y="A549")

ggsave("plots/scatter_facs_celllines.pdf",width = 6, height = 4)
