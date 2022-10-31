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

#%% Load data
codonocc = pd.read_csv("results/all_CodPairs_ratios_UPDOWNproteins.csv",index_col=0)

# Remove NAs
codonocc.replace([np.inf, -np.inf], np.nan, inplace=True)
codonocc.dropna(axis=1,how="all", inplace=True)
data = codonocc.replace(np.nan, 0)
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
codclusters = {"cluster1":["Uterus","Smoothmuscle","Breast","Vagina","Tonsil","Colon",
                       "Thyroid","Muscle...Skeletal","Artery","Fat","Heart","Esophagus",
                       "Ovary","Liver","Lymphnode","Fallopiantube","Rectum","Kidney"],
               "cluster2":["Lung","Appendices","Placenta","Brain","Adrenal","Testis",
                       "Stomach","Skin","Nerve...Tibial","Gallbladder","Urinarybladder",
                       "Pituitary","Smallintestine","Duodenum","Prostate","Spleen",
                       "Salivarygland","Pancreas"]}
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(codclusters))/(len(codclusters)-1))
handles = [mlines.Line2D([], [], color=c, alpha=0.5) for c in colors]

for label, color in zip(codclusters.keys(),colors):
    idx = [s in codclusters[label] for s in principalDf.index]
    ax.scatter(principalDf.loc[idx, 'PCA1'], principalDf.loc[idx, 'PCA2'], c = color, marker = 'o', s = 50, alpha=0.5)
    for lab in principalDf.index[idx]:
        ax.annotate(lab,[principalDf.loc[lab, 'PCA1'],principalDf.loc[lab, 'PCA2']],fontsize=8)
leg = fig.legend(handles,codclusters.keys(),bbox_to_anchor=(0.37,0.2,0.5,0.5))
fig.add_artist(leg)
ax.grid()

#%%
features.to_csv("results/PCA_ratios_CodPairs.csv")