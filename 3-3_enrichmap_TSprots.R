
# Load tissue-specificities
TSprots = read.delim("results/TSall_PTR-UP-FCavg2-FCglob1.tsv",sep="\t")

### Generate dataframe file for Cytoscape ###
outUP = data.frame(row.names = seq(length(tissues)))
outUP$genesetID = tissues
outUP$description = tissues
outUP$pvalue = 0
outUP$FDR = 0
outUP$Phenotype = 1
outUP$genelist = sapply(tissues, function(x) paste(up_list[[x]],collapse = ","))
write.table(outUP,"results/UPall_PTR-UP-FCavg2-FCglob1.txt",sep="\t",quote=F,row.names = F)
outDOWN = data.frame(row.names = seq(length(tissues)))
outDOWN$genesetID = tissues
outDOWN$description = tissues
outDOWN$pvalue = 0
outDOWN$FDR = 0
outDOWN$Phenotype = -1
outDOWN$genelist = sapply(tissues, function(x) paste(down_list[[x]],collapse = ","))
write.table(outDOWN,"results/DOWNall_PTR-UP-FCavg2-FCglob1.txt",sep="\t",quote=F,row.names = F)
