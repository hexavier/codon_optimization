#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import seaborn as sns
def norm_len(codus):
    codnorm = pd.Series(index = codus.index, dtype="float64")
    total = codus.sum()
    for cod in codus.index:
        codnorm.loc[cod] = codus.loc[cod]/total
    return codnorm

#%% Create dataframe with all features and targets
# Upload data
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
codons = list(GENETIC_CODE.keys())
codons.remove("TAA"); codons.remove("TGA"); codons.remove("TAG")

codus = pd.read_csv("data/human_CU_ENS_09062020.tsv", sep="\t", index_col=0)
TSproteins = pd.read_csv("results/TSall_PTR-UP-FCavg2-FCglob1.tsv", sep="\t", index_col=0)
HLproteins = pd.read_csv("data/ProteinHalfLife_merged.csv", index_col=0).mean(axis=1)

#%% For each tissue, compute ratios

# Create results object
tissues = list(set(TSproteins.Tissue))
ratios = pd.DataFrame(index = tissues, columns = codons, dtype="float64")
prct = 40

for t in tissues:
    TStemp = TSproteins.loc[TSproteins.Tissue==t,:]
    genes = [s for s in codus.index if (s in TStemp.index) and (s in HLproteins.index)]
    codusnorm = codus.loc[genes,codons].apply(norm_len, axis=1)
    indf = pd.concat([codusnorm,TStemp.loc[genes,["updown"]],HLproteins[genes]], axis=1)
    # Among the genes that are tissue specific, take those with high PTR and low HL
    # versus low PTR and high HL
    UP = indf.loc[indf.updown=="UP",]
    UP = UP.loc[UP[0]<np.percentile(UP[0],prct),codons]
    DOWN = indf.loc[indf.updown=="DOWN",]
    DOWN = DOWN.loc[DOWN[0]>np.percentile(DOWN[0],100-prct),codons]
    print("%s: %i proteins UP and %i proteins DOWN" % (t,UP.shape[0],DOWN.shape[0]))
    ratios.loc[t,codons] = np.log2(UP.mean()/DOWN.mean())
    
#%% Plot
g = sns.clustermap(ratios, center=0, cmap="vlag", method="complete",figsize=(16,10))

