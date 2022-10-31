## Thresholds
nfold_glob = 1
nfold_avg = 2

## Load PTR data
# GTEX data
tissuePTR = read.csv("data/PTR_medians_GTEx_Jiang2020.csv", row.names = 1)

# Define TS proteins
TSevents = c()
for (tissue in colnames(tissuePTR)){
  # Compute metrics to define tissue-specificity
  maxothers = apply(tissuePTR,1,function(x) max(x[colnames(tissuePTR)!=tissue], na.rm = T))
  minothers = apply(tissuePTR,1,function(x) min(x[colnames(tissuePTR)!=tissue], na.rm = T))
  Nothers = rowSums(!is.na(tissuePTR[,colnames(tissuePTR)!=tissue]))
  avgothers = apply(tissuePTR,1,function(x) mean(x[colnames(tissuePTR)!=tissue], na.rm = T))
  upFCglob = tissuePTR[,tissue] - maxothers
  downFCglob = tissuePTR[,tissue] - minothers
  FCavg = tissuePTR[,tissue] - avgothers
  
  # Find TS events. Note that PTR data is log-10 transformed
  isUP = (upFCglob>= log10(nfold_glob))&(FCavg>=log10(nfold_avg))&(Nothers>=3)
  isUP[is.na(isUP)] <- F
  
  isDOWN = (downFCglob<= -log10(nfold_glob))&(FCavg<=-log10(nfold_avg))&(Nothers>=3)
  isDOWN[is.na(isDOWN)] <- F
  
  # Create temp dataset
  tempevents = data.frame(row.names = rownames(tissuePTR[isUP|isDOWN,]))
  tempevents$gene = rownames(tissuePTR[isUP|isDOWN,])
  tempevents$Tissue = tissue
  tempevents$PTRtissue = tissuePTR[isUP|isDOWN,tissue]
  tempevents$avgothers = avgothers[isUP|isDOWN]
  tempevents$FCavg = FCavg[isUP|isDOWN]
  tempevents$maxothers = maxothers[isUP|isDOWN]
  tempevents$upFCglob = upFCglob[isUP|isDOWN]
  tempevents$minothers = minothers[isUP|isDOWN]
  tempevents$downFCglob = downFCglob[isUP|isDOWN]
  tempevents$Nothers = Nothers[isUP|isDOWN]
  tempevents$updown = "UP"
  tempevents[rownames(tissuePTR[isDOWN,]),"updown"] = "DOWN"

  # Merge
  TSevents = rbind(TSevents,tempevents,make.row.names=T)
}

# Get consensus list of proteins
TSmergedprots = c()
for (t in unique(TSevents$Tissue)){
  tisprots = TSevents[TSevents$Tissue==t,]
  TSmergedprots = rbind(TSmergedprots,
                        unique(tisprots[,c("gene","Tissue","updown")]))
}


# Save
write.table(TSmergedprots,sprintf("results/TSgtex_PTR-UP-FCavg%s-FCglob%s.tsv",nfold_avg,nfold_glob),sep="\t",row.names = F)
