## Thresholds
nfold_glob = 0
nfold_avg = 2

## Load PTR data
# GTEX data
gtex = read.csv("data/PTR_medians_GTEx_Jiang2020.csv", row.names = 1)
# Eraslan(2019) data
eraslan = read.csv("data/PTR_Eraslan2019.tsv", sep="\t", row.names = 2)
eraslan = eraslan[,4:ncol(eraslan)]

# Define TS proteins
TSevents = c()
for (d in c("gtex","eraslan")){
  tissuePTR = get(d)
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
    tempevents$dataset = d
    
    # Merge
    TSevents = rbind(TSevents,tempevents,make.row.names=T)
  }
}

# Get consensus list of proteins
TSmergedprots = c()
for (t in unique(TSevents$Tissue)){
  tisprots = TSevents[TSevents$Tissue==t,]
  # Identify unconsistent changes
  up_eraslan = tisprots[(tisprots$updown=="UP")&(tisprots$dataset=="eraslan"),"gene"]
  down_eraslan = tisprots[(tisprots$updown=="DOWN")&(tisprots$dataset=="eraslan"),"gene"]
  up_gtex = tisprots[(tisprots$updown=="UP")&(tisprots$dataset=="gtex"),"gene"]
  down_gtex = tisprots[(tisprots$updown=="DOWN")&(tisprots$dataset=="gtex"),"gene"]
  inconsistent = c(up_eraslan[up_eraslan %in% down_gtex],up_gtex[up_gtex %in% down_eraslan])
  # Keep union of both datasets
  TSmergedprots = rbind(TSmergedprots,
                        unique(tisprots[!(tisprots$gene %in% inconsistent),c("gene","Tissue","updown")]))
}


# Save
write.table(TSmergedprots,sprintf("results/TSall_PTR-UP-FCavg%s-FCglob%s.tsv",nfold_avg,nfold_glob),sep="\t",row.names = F)

### Check overlap ###
matched = colnames(eraslan)[colnames(eraslan) %in% colnames(gtex)]
overlap_df = data.frame(row.names = matched)
for (t in matched){
  tissue_exp = TSevents[TSevents$Tissue==t,]
  tissue_gtex_up = tissue_exp[(tissue_exp$dataset=="gtex")&(tissue_exp$updown=="UP"),"gene"]
  tissue_eraslan_up = tissue_exp[(tissue_exp$dataset=="eraslan")&(tissue_exp$updown=="UP"),"gene"]
  tissue_gtex_down = tissue_exp[(tissue_exp$dataset=="gtex")&(tissue_exp$updown=="DOWN"),"gene"]
  tissue_eraslan_down = tissue_exp[(tissue_exp$dataset=="eraslan")&(tissue_exp$updown=="DOWN"),"gene"]
  overlap_df[t,"agreed"] = sum(tissue_gtex_up %in% tissue_eraslan_up) + sum(tissue_gtex_down %in% tissue_eraslan_down)
  overlap_df[t,"nonagreed"] = sum(tissue_gtex_down %in% tissue_eraslan_up) + sum(tissue_gtex_up %in% tissue_eraslan_down)
  overlap_df[t,"total"] = length(unique(tissue_exp$gene))
  overlap_df[t,"binomial_test"] = binom.test(overlap_df[t,"agreed"],overlap_df[t,"agreed"]+overlap_df[t,"nonagreed"],alternative="greater")$p.value
}
write.csv(overlap_df,sprintf("results/overlap_datasets_TSproteins-FCavg%s-FCglob%s.csv",nfold_avg,nfold_glob))
