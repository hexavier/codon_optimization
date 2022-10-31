## Load data
cocoputs = read.csv("data/o586358-Human_CDS_Bicod_25012021.tsv",sep="\t") # from https://hive.biochemistry.gwu.edu/review/codon2
mapping = read.csv("data/MANE.GRCh38.v0.91.summary.txt",sep="\t")

## Map RefSeq IDs to ENSembl
# Remove version from CoCoPUTs ID
cocoputs$Protein.ID = as.character(sapply(as.character(cocoputs$Protein.ID), function(x) strsplit(x,"\\.")[[1]][1]))

# Remove version from mapping IDs
mapping$Ensembl_Gene = as.character(sapply(as.character(mapping$Ensembl_Gene), function(x) strsplit(x,"\\.")[[1]][1]))
mapping$RefSeq_prot = as.character(sapply(as.character(mapping$RefSeq_prot), function(x) strsplit(x,"\\.")[[1]][1]))
# Get average CU of different proteins from same gene
cocoputs_ens = sapply(as.character(cocoputs$Protein.ID), 
                      function(x) if(x %in% mapping$RefSeq_prot){mapping[mapping$RefSeq_prot %in% x,"Ensembl_Gene"]}else{NA})
codus = sapply(unique(cocoputs_ens)[!is.na(unique(cocoputs_ens))], 
               function(x) colMeans(cocoputs[cocoputs_ens %in% x,9:4110]))

## Save
write.table(t(codus), "data/human_CodPair_ENS_25012021.tsv",sep="\t")
