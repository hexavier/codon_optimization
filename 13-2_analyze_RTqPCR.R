library(ggplot2)
library(ggpubr)

# Load data
ct_data = read.delim("data/qPCR_summary.tsv", sep="\t")

# Write full names
ct_data$construct = factor(as.character(ct_data$construct), levels=c("mCLeGK1","mCLeGK2","mCKeGL1","mCKeGL2"))
ct_data$cell = factor(as.character(ct_data$cell), levels=c("A549","HEK293"))

# Plot
ggplot(ct_data, aes(x=construct,y=dCt,fill=cell, group=cell)) +
  facet_wrap(. ~ construct, ncol = 4, scales = "free") +
  geom_point(shape=21, color="black", size=3, position=position_dodge(width=0.6)) +
  scale_fill_manual(values=c("#304d63","#8fb9aa")) +
  expand_limits(y = 0) +
  theme_classic() +
  labs(y = "2^-dCt")

ggsave("plots/rtqpcr.pdf",width = 8, height = 2.5)
