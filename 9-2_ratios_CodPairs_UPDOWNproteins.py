#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
def norm_len(codus):
    codnorm = pd.Series(index = codus.index, dtype="float64")
    total = codus.sum()
    for cod in codus.index:
        codnorm.loc[cod] = codus.loc[cod]/total
    return codnorm

#%% Create dataframe with all features and targets
# Upload data
codus = pd.read_csv("data/human_CodPair_ENS_25012021.tsv", sep="\t", index_col=0)
TSproteins = pd.read_csv("results/TSall_PTR-UP-FCavg2-FCglob1.tsv", sep="\t", index_col=0)

#%% For each tropism, create a RF model that distinguishes one tropism from others based on SDA

# Create results object
tissues = list(set(TSproteins.Tissue))
codons = list(codus.columns)[5:]
ratios = pd.DataFrame(index = tissues, columns = codons)

for t in tissues:
    TStemp = TSproteins.loc[TSproteins.Tissue==t,:]
    genes = [s for s in codus.index if s in TStemp.index]
    codusnorm = codus.loc[genes,codons].apply(norm_len, axis=1)
    indf = pd.concat([codusnorm,TStemp.loc[genes,["updown"]]], axis=1)
    means = indf.groupby("updown").mean()
    ratios.loc[t,codons] = np.log2(means.loc["UP",codons]/means.loc["DOWN",codons])

#%% Save
ratios.to_csv("results/all_CodPairs_ratios_UPDOWNproteins.csv")

#%% Check numbers of UP and DOWN codon pairs

# Create results object
tissues = list(set(TSproteins.Tissue))
codons = list(codus.columns)[5:]
UPcounts = pd.DataFrame(index = tissues, columns = codons)
DOWNcounts = pd.DataFrame(index = tissues, columns = codons)

for t in tissues:
    TStemp = TSproteins.loc[TSproteins.Tissue==t,:]
    genes = [s for s in codus.index if s in TStemp.index]
    codustemp = codus.loc[genes,codons]
    indf = pd.concat([codustemp,TStemp.loc[genes,["updown"]]], axis=1)
    sums = indf.groupby("updown").sum()
    UPcounts.loc[t,codons] = sums.loc["UP",codons]
    DOWNcounts.loc[t,codons] = sums.loc["DOWN",codons]

UPcounts.to_csv("results/all_CodPairs_counts_UPproteins.csv")
DOWNcounts.to_csv("results/all_CodPairs_counts_DOWNproteins.csv")
