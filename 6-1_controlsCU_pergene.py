# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

def get_codons(seq,frame):
    cutseq = seq[frame:]
    seqcodons = (cutseq[n:n+3] for n in range(0,len(cutseq),3) if len(cutseq)>=(n+3))
    return seqcodons
def get_dinucl(seq):
    seqdinucl = (seq[n:n+2] for n in range(0,len(seq),2) if len(seq)>=(n+2))
    return seqdinucl
def openfasta(filename):
    # Open the file
    fileHandle = open(filename)
    # Read all the lines in the file
    seqs = fileHandle.read().split(">") # Split file whenever a > is detected
    # Close the file
    fileHandle.close()
    # Loop over the lines
    seqdict = {}
    for seq in seqs:
        if seq: # Skip empty sequences
            splitted = seq.split("\n") # Split different lines of one same sequence
            # the first line is the name of the sequence, the rest is the sequence itself
            seqdict[splitted[0]] = "".join(splitted[1:])
    return seqdict

GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

DINUCL = ["AA","TT","GG","CC","AT","AG","AC","TA","TC","TG","GA","GC","GT","CA","CT","CG"]

#%% Load data
# Fasta sequences parsing
ccds_seq = openfasta("data/CCDS_nucleotide.current.fna")
# CCDSs
ccds = pd.read_csv("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt", sep="\t", index_col=4)


#%% Create output table
# Find gene names
ccds_id = {g.split("|")[0]:g for g in ccds_seq.keys()}

# Initialize structures
frame1 = pd.DataFrame(columns = GENETIC_CODE.keys(), index = ccds.index, dtype=float)
frame2 = pd.DataFrame(columns = GENETIC_CODE.keys(), index = ccds.index, dtype=float)
dinucl = pd.DataFrame(columns = DINUCL, index = ccds.index, dtype=float)

# Count codons
for p in frame1.index:
    if p in ccds_id.keys():
        # Find sequence
        seq = ccds_seq[ccds_id[p]]
        
        # Count codons frame +1
        seqcodons = get_codons(seq,1)
        frame1.loc[p,:] = 0
        for c in seqcodons:
            if c in frame1.columns:
                frame1.loc[p,c] += 1
                
        # Count codons frame +2
        seqcodons = get_codons(seq,2)
        frame2.loc[p,:] = 0
        for c in seqcodons:
            if c in frame2.columns:
                frame2.loc[p,c] += 1
                
        # Count dinucleotides
        seqcodons = get_dinucl(seq)
        dinucl.loc[p,:] = 0
        for c in seqcodons:
            if c in dinucl.columns:
                dinucl.loc[p,c] += 1
    else:
        frame1.loc[p,:] = np.nan
        frame2.loc[p,:] = np.nan
        dinucl.loc[p,:] = np.nan

frame1.dropna(inplace=True)
frame2.dropna(inplace=True)
dinucl.dropna(inplace=True)
#%% Get average per gene
# Mapping files
prot2gene = pd.read_csv("data/MANE.GRCh38.v0.91.summary.txt",sep="\t")
prot2gene["Ensembl_Gene"] = [s.split(".")[0] for s in prot2gene["Ensembl_Gene"]]
prot2gene["RefSeq_prot"] = [s.split(".")[0] for s in prot2gene["RefSeq_prot"]]
ccds2prot = pd.read_csv("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS2Sequence.current.txt",sep="\t",index_col="#ccds")
ccds2prot = ccds2prot.loc[ccds2prot.source=="NCBI"]
# Identify genes
genes = pd.DataFrame(index = frame1.index, columns=["gene"])
for p in frame1.index:
    prots = ccds2prot.loc[ccds2prot.index==p,"protein_ID"].unique()
    for n in prots:
        prot = n.split(".")[0]
        if prot in prot2gene.RefSeq_prot.values:
            genes.loc[p,"gene"] = prot2gene.loc[prot2gene.RefSeq_prot==prot,"Ensembl_Gene"].values[0]

frame1 = pd.concat([frame1,genes],axis=1); frame1.dropna(inplace=True)
frame2 = pd.concat([frame2,genes],axis=1); frame2.dropna(inplace=True)
dinucl = pd.concat([dinucl,genes],axis=1); dinucl.dropna(inplace=True)
# Get average per gene
frame1_gene = frame1.groupby("gene").mean()
frame2_gene = frame2.groupby("gene").mean()
dinucl_gene = dinucl.groupby("gene").mean()

#%% Save output
frame1_gene.to_csv("data/human_CU_ENS_frame1.tsv", sep="\t")
frame2_gene.to_csv("data/human_CU_ENS_frame2.tsv", sep="\t")
dinucl_gene.to_csv("data/human_CU_ENS_dinucl.tsv", sep="\t")
