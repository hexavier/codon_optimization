# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 09:11:58 2020

@author: user
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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

def relative_codons(codus):
    codnorm = pd.Series(index = codus.index, dtype="float64")
    AAs = list(set(GENETIC_CODE.values())); AAs.remove("_")
    for aa in AAs:
        cod_aa = [c for c in GENETIC_CODE.keys() if GENETIC_CODE[c]==aa]
        total_aa = codus.loc[cod_aa].sum()
        if total_aa > 0.0:
            for cod in cod_aa:
                codnorm.loc[cod] = codus.loc[cod]/total_aa
        else:
            for cod in cod_aa:
                codnorm.loc[cod] = 1.0/len(cod_aa)
    return codnorm

#%% Open fasta
proteins = openfasta("results/optimized_fluorescent_proteins.fa")

#%% Compute RCUs
codoncount = pd.DataFrame(index=proteins.keys(), columns = GENETIC_CODE.keys())
for p in codoncount.index:
    seq = proteins[p]
    seqcodons = [seq[n:n+3] for n in range(0,len(seq),3)]
    codoncount.loc[p] = [seqcodons.count(cod) for cod in codoncount.columns]
rcu = codoncount.apply(relative_codons,1)

#%% Plot
diff = np.subtract(rcu.iloc[[0,2],:],rcu.iloc[[1,3],:])

diff.T.plot.bar(figsize=(10,5))

diff.T.plot.scatter(x=0,y=1)
for label in diff.columns:
    plt.annotate(label,[diff[label].iloc[0], diff[label].iloc[1]], fontsize=7)

correlations = diff.T.corr(method="spearman")

#%% Plot heatmap of composition
codusnormratios = pd.read_csv("results/all_codusnorm_ratios_UPDOWNproteins.csv", index_col = 0)
for p in proteins.keys():
    for tissue in ["Kidney","Lung"]:
        seq = proteins[p]
        seqcodons = [seq[n:n+3] for n in range(0,len(seq),3)]
        codon_weights = np.array([codusnormratios.loc[tissue,c] for c in seqcodons if c in codusnormratios.columns])
        codcomposition = pd.DataFrame(codon_weights.reshape(1,len(codon_weights)),
                                      index=[p],columns=range(len(codon_weights)))
        # Plot
        fig, ax = plt.subplots(figsize=(10,0.5))
        sns.heatmap(codcomposition, vmin=-1, vmax=1, cmap="bwr", xticklabels=False, yticklabels=False)
        fig.savefig(str("plots/%s_to%s.pdf" % (p,tissue)))
#%% Write rcu
rcu.to_csv("results/RCUs_optimized_proteins.csv")