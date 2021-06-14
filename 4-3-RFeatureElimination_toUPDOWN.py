#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

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

#%% For each tissue, run RF to distinguish UP vs DOWN

tissues = list(set(TSproteins.Tissue))
# Create results object
AUCmeans = pd.DataFrame(index = tissues, columns = np.arange(len(codons)))
AUCstd = pd.DataFrame(index = tissues, columns = np.arange(len(codons)))
RFEmeans = pd.DataFrame(index = tissues, columns = list(codons))
RFEstd = pd.DataFrame(index = tissues, columns = list(codons))

for t in tissues:
    TStemp = TSproteins.loc[TSproteins.Tissue==t,:]
    genes = [s for s in codus.index if s in TStemp.index]
    codusnorm = codus.loc[genes,codons].apply(norm_len, axis=1)
    indf = pd.concat([codusnorm,TStemp.loc[genes,["updown"]]], axis=1)
    X = indf.iloc[:,0:-1]  # Features
    y = indf.iloc[:,-1]=="UP"  # Labels: 1 for UP, 0 for DOWN
    # Check number of target groups
    print("%s: There are %i proteins UP, and %i proteins DOWN" % (t,sum(y),(len(y)-sum(y))))
    # To make groups comparable, separate equal samples of target and non-target
    roc_aucs = []; ranks = [];
    n_limiting = min(sum(y),(len(y)-sum(y))) # Check whether limiting group is target or non-target
    for n in range(100):
        # Separe targets from non-targets and sample from them
        targets = indf.loc[y,:].sample(n=n_limiting)
        nontargets = indf.loc[~y,:].sample(n=n_limiting)
        tempdf = targets.append(nontargets)
        # Define features and targets
        Xtemp = tempdf.iloc[:,0:-1]
        ytemp = tempdf.iloc[:,-1]=="UP"
        # Create a Random Forest Classifier
        clf=RandomForestClassifier(n_estimators=100)
        # Initialize RFECV object
        feature_selector = RFECV(clf, cv = 5, step = 1, scoring = "roc_auc")
        # Fit RFECV
        feature_selector.fit(Xtemp, ytemp)
        # Get selected features
        feature_names = Xtemp.columns
        roc_aucs.append(feature_selector.grid_scores_)
        ranks.append(feature_selector.ranking_)
    
    # Record results
    AUCmeans.loc[t,:] = np.mean(roc_aucs,axis=0)
    AUCstd.loc[t,:] = np.std(roc_aucs,axis=0)
    RFEmeans.loc[t,:] = np.mean(ranks,axis=0)
    RFEstd.loc[t,:] = np.std(ranks,axis=0)
    
#%% #ROC curve
    fig = plt.figure()
    mean_auc = np.mean(roc_aucs,axis=0)
    std_auc = np.std(roc_aucs,axis=0)
    aucs_upper = np.minimum(mean_auc + std_auc, 1)
    aucs_lower = np.maximum(mean_auc - std_auc, 0)
    plt.plot(np.arange(1,len(codons)+1), mean_auc, marker='.', color="b", lw = 2)
    plt.fill_between(np.arange(1,len(codons)+1), aucs_lower, aucs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')
    # axis labels
    plt.title(str("RFE of %s" % t))
    plt.xlabel('Number of Features')
    plt.ylabel('AUC')
    # show the legend
    plt.legend(loc="lower right")
    # show the plot
    plt.show()
    fig.savefig(str("plots/all_RFE_codusnorm_to_UPDOWNproteins/RFE_%s.pdf" % t), bbox_inches='tight')

#%% #ROC curve
cmap = plt.cm.get_cmap('jet') # set color map
colors = cmap(np.arange(len(tissues))/(len(tissues)-1))
handles = [mlines.Line2D([], [], marker='o', color=c, alpha=0.5) for c in colors]

fig = plt.figure()
for t, color in zip(tissues,colors):
    plt.plot(np.arange(1,len(codons)+1), AUCmeans.loc[t,:], marker='.', color=color, lw = 2,
             label= t)

# axis labels
plt.title("RFE of all tissues")
plt.xlabel('Number of Features')
plt.ylabel('AUC')
# show the legend
fig.legend(handles,tissues,loc='center left', bbox_to_anchor=(0.95, 0.5))
# show the plot
plt.show()
fig.savefig("plots/all_RFE_codusnorm_to_UPDOWNproteins/RFE.pdf", bbox_inches='tight')
    
#%% Save
RFEmeans.to_csv("results/mean_all_RFE_codusnorm_to_UPDOWNproteins.csv")
RFEstd.to_csv("results/std_all_RFE_codusnorm_to_UPDOWNproteins.csv")
