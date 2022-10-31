library(ggplot2)
library(ggpubr)
library(patchwork)
library(reshape)
library(ggrepel)

# Load data
raw_data = read.csv("data/2022MC004_GFP_mcherry_001to024_table-results_NoiseIntegration.csv")

# Take average 
avg_data = cast(raw_data, Sample.name ~ Protein, mean, value="Normalized.Area", na.rm=T)
avg_data$ratio = log2(avg_data$GFP / avg_data$mCherry)
avg_data$plasmid = sapply(avg_data$Sample.name,function(x) strsplit(x,"-|_")[[1]][2])
avg_data$cell = sapply(avg_data$Sample.name,function(x) strsplit(x,"-|_")[[1]][1])

# Write full names
avg_data$plasmid = sapply(avg_data$plasmid, function(x) if(x=="74"){"mCherry_L-eGFP_K-1"}else 
  if(x=="75"){"mCherry_K-eGFP_L-1"}else 
    if(x=="76"){"mCherry_L-eGFP_K-2"} else 
      if(x=="77"){"mCherry_K-eGFP_L-2"})
avg_data$plasmid = factor(as.character(avg_data$plasmid), levels=c("mCherry_L-eGFP_K-1","mCherry_L-eGFP_K-2","mCherry_K-eGFP_L-1","mCherry_K-eGFP_L-2"))
avg_data$cell = factor(as.character(avg_data$cell), levels=c("A549","293T","lung","renal"))

# Separate cell lines and primary
celllines = avg_data[avg_data$cell %in% c("293T", "A549"),]
primary = avg_data[avg_data$cell %in% c("lung", "renal"),]

# Plot
p1 <- ggplot(celllines, aes(x=plasmid,y=ratio,fill=cell, group=cell)) +
  geom_hline(yintercept=0) +
  geom_point(shape=21, color="black", size=3, position=position_dodge(width=0.6)) +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  theme_classic() +
  labs(y = "log2(eGFP/mCherry)")

p2 <- ggplot(primary, aes(x=plasmid,y=ratio,fill=cell, group=cell)) +
  geom_hline(yintercept=0) +
  geom_point(shape=21, color="black", size=3, position=position_dodge(width=0.6)) +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  theme_classic() +
  labs(y = "log2(eGFP/mCherry)")

p1 + p2 + plot_layout(widths = c(2, 1))
ggsave("plots/proteomics.pdf",width = 10, height = 2.5)

#### Check absolute abundances ####
p1 <- ggplot(celllines, aes(x=GFP,y=mCherry,fill=cell, label=plasmid)) +
  geom_point(shape=21, color="black", size=3) +
  geom_text_repel(size=3) +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  scale_y_log10() + scale_x_log10() +
  theme_classic()
p2 <- ggplot(primary, aes(x=GFP,y=mCherry,fill=cell, label=plasmid)) +
  geom_point(shape=21, color="black", size=3) +
  geom_text_repel(size=3) +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  scale_y_log10() + scale_x_log10() +
  theme_classic()
p1|p2
ggsave("plots/proteomics_scatter.pdf",width = 12, height = 5)
