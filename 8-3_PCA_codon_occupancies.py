#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 16:09:57 2019

@author: xhernandez
"""

# Load modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.lines as mlines
from combat.pycombat import pycombat

#%% Load data
codonocc = pd.read_csv("data/codon_occupancies.csv",index_col=0)
mapping = pd.read_table("data/riboseq_datasets.tsv", index_col=0)
# Remove stop codons
codonocc = codonocc.loc[:,[c not in ["TAA","TAG","TGA"] for c in codonocc.columns]]

#%% Correct batches with ComBat
notcorr = codonocc.dropna().T
mapbatch = dict(zip(list(set(mapping.study)),[0,1,2,3,4]))
batch = [mapbatch[mapping.loc[s,"study"]] for s in notcorr.columns]
data = pycombat(notcorr,batch).T

#%% PCA
pca = PCA(n_components=2)
x = StandardScaler().fit_transform(data)
principalComponents = pca.fit_transform(x)
expl_var = pca.explained_variance_ratio_
features = pd.DataFrame(pca.components_.transpose(), columns=["PCA1","PCA2"], index = data.columns)

## Plot
principalDf = pd.DataFrame(principalComponents,columns=["PCA1","PCA2"],index=data.index)
 
# Plot pca
fig = plt.figure(figsize = (10,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel(str('Dim 1 (%1.2f%%)' % (expl_var[0]*100)), fontsize = 15)
ax.set_ylabel(str('Dim 2 (%1.2f%%)' % (expl_var[1]*100)), fontsize = 15)
ax.set_title('Dimension Reduction', fontsize = 20)

# Color based on tissue, shape based on study
tissues = list(set(mapping.tissue))
datasets = dict(zip(list(set(mapping.study)),["o","v","d","*","s"]))
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(tissues))/(len(tissues)-1))
handles1 = [mlines.Line2D([], [], color=c, alpha=0.5) for c in colors]
handles2 = [mlines.Line2D([], [], marker=datasets[s], alpha=0.5) for s in datasets.keys()]

for label, color in zip(tissues,colors):
    for d,shape in datasets.items():
        idx = list([np.logical_and(mapping.loc[s,"tissue"]==label,mapping.loc[s,"study"]==d) for s in principalDf.index])
        ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, marker = shape, s = 50, alpha=0.5)
        for lab in principalDf.index[idx]:
           ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
leg1 = fig.legend(handles1,tissues,bbox_to_anchor=(0.37,0.2,0.5,0.5))
leg2 = fig.legend(handles2,datasets.keys(),bbox_to_anchor=(0.37,0.35,0.5,0.5))
fig.add_artist(leg1)
ax.grid()

#%%
features.to_csv("results/PCA_ComBat_codon_occupancies.csv")