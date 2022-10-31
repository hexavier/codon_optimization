#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import auc
import matplotlib.pyplot as plt
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
TSproteins = pd.read_csv("results/TSgtex_PTR-UP-FCavg2-FCglob1.tsv", sep="\t", index_col=0)

#%% For each tissue, run RF to distinguish UP vs DOWN

# Create results object
colnames = ["roc_auc","accuracy"]; colnames.extend(codons)
tissues = list(set(TSproteins.Tissue))
RFmeans = pd.DataFrame(index = tissues, columns = colnames)
RFstd = pd.DataFrame(index = tissues, columns = colnames)

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
    accuracies = []; roc_aucs = []; prec_aucs = []; 
    features = pd.DataFrame(index=indf.iloc[:,0:-1].columns)
    n_limiting = min(sum(y),(len(y)-sum(y))) # Check whether limiting group is target or non-target
    # Initiate structures for ROC curve
    tprs = []; aucs = [] ; precisions = []
    linear_xaxis = np.linspace(0, 1, 100)
    for n in range(100):
        # Separe targets from non-targets and sample from them
        targets = indf.loc[y,:].sample(n=n_limiting)
        nontargets = indf.loc[~y,:].sample(n=n_limiting)
        tempdf = pd.concat([targets,nontargets])
        # Define features and targets
        Xtemp = tempdf.iloc[:,0:-1]
        ytemp = tempdf.iloc[:,-1]=="UP"
        # Create a Random Forest Classifier
        clf=RandomForestClassifier(n_estimators=100)
        # Cross validation with 5 splits
        cv = StratifiedKFold(n_splits=5)
        acctemp = []; roc_auctemp = []; feattemp = []; tprtemp = []; prectemp = []; prec_auctemp = []
        for i, (train, test) in enumerate(cv.split(Xtemp, ytemp)):
            clf.fit(Xtemp.iloc[train,:], ytemp.iloc[train])
            # ROC curve and Precision-Recall curve
            # predict probabilities
            clf_probs = clf.predict_proba(Xtemp.iloc[test,:])
            # keep probabilities for the positive outcome only
            clf_probs = clf_probs[:, 1]
            fpr, tpr, _ = roc_curve(ytemp.iloc[test], clf_probs)
            precision, recall, _ = precision_recall_curve(ytemp.iloc[test], clf_probs)
            # Interpolate tpr values to match a fixed x axis
            interp_tpr = np.interp(linear_xaxis, fpr, tpr) 
            interp_tpr[0] = 0.0
            interp_precision = np.interp(linear_xaxis, precision, recall) 
            interp_precision[0] = 1.0
            # Save iteration
            tprtemp.append(interp_tpr)
            prectemp.append(interp_precision)
            roc_auctemp.append(roc_auc_score(ytemp.iloc[test], clf_probs))
            prec_auctemp.append(auc(recall,precision))
            # Accuracy
            y_pred=clf.predict(Xtemp.iloc[test,:])
            acctemp.append(accuracy_score(ytemp.iloc[test],y_pred))
            # Features
            feattemp.append(clf.feature_importances_)
        
        # Model Accuracy, how often is the classifier correct?
        roc_aucs.append(np.mean(roc_auctemp))
        prec_aucs.append(np.mean(prec_auctemp))
        tprs.append(pd.DataFrame(tprtemp).mean(axis = 0))
        precisions.append(pd.DataFrame(prectemp).mean(axis = 0))
        accuracies.append(np.mean(acctemp))
        # Analyze features
        features.loc[:,n] = pd.DataFrame(feattemp).mean(axis=0).values
    
    # Record results
    trop_means = [np.mean(roc_aucs), np.mean(accuracies)]; trop_means.extend(features.mean(axis=1))
    RFmeans.loc[t,:] = trop_means
    trop_std = [np.std(roc_aucs), np.std(accuracies)]; trop_std.extend(features.std(axis=1))
    RFstd.loc[t,:] = trop_std
    
#%% #Creating a bar plot of feature importances
    fig = plt.figure(figsize=(10,12))
    features = features.reindex(features.mean(axis=1).sort_values(ascending=False).index, axis=0)
    sns.barplot(x=features.mean(axis=1), y=features.index, ci=features.std(axis=1))
    plt.errorbar(x= features.mean(axis=1), y= features.index, fmt='none', xerr=features.std(axis=1), ecolor="k")
    # Add labels to your graph
    plt.xlabel('Feature Importance Score')
    plt.ylabel('Features')
    plt.title(str("Visualizing Important Features in %s" % t))
    plt.show()
    fig.savefig(str("plots/gtex_RandomForest_codusnorm_to_UPDOWNproteins/features_proteins_%s.pdf" % t), bbox_inches='tight')
#%% #ROC curve
    fig = plt.figure()
    mean_tpr = np.mean(tprs, axis=0); std_tpr = np.std(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = np.mean(roc_aucs); std_auc = np.std(roc_aucs)
    
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
    plt.plot(linear_xaxis, mean_tpr, marker='.', color="b", lw = 2,
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc))
    plt.fill_between(linear_xaxis, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')
    # axis labels
    plt.title(str("ROC curve of %s" % t))
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    # show the legend
    plt.legend(loc="lower right")
    # show the plot
    plt.show()
    fig.savefig(str("plots/gtex_RandomForest_codusnorm_to_UPDOWNproteins/roc_proteins_%s.pdf" % t), bbox_inches='tight')
    
#%% #Precision-Recall Curve
    fig = plt.figure()
    mean_prec = np.mean(precisions, axis=0); std_prec = np.std(precisions, axis=0)
    mean_prec[-1] = 0.0
    mean_auc = np.mean(prec_aucs); std_auc = np.std(prec_aucs)
    
    precs_upper = np.minimum(mean_prec + std_prec, 1)
    precs_lower = np.maximum(mean_prec - std_prec, 0)
    
    chance = len(y[y==1]) / len(y)
    plt.plot([0, 1], [chance, chance], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)
    plt.plot(linear_xaxis, mean_prec, marker='.', color="b", lw = 2,
             label=r'Mean Precision-Recall (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc))
    plt.fill_between(linear_xaxis, precs_lower, precs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')
    # axis labels
    plt.title(str("Precision-Recall Curve curve of %s" % t))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    # show the legend
    plt.legend(loc="lower left")
    # show the plot
    plt.show()
    fig.savefig(str("plots/gtex_RandomForest_codusnorm_to_UPDOWNproteins/prec-recall_proteins_%s.pdf" % t), bbox_inches='tight')

#%% Save
RFmeans.to_csv("results/mean_ptr-gtex_randomforest_codusnorm_to_UPDOWNproteins.csv")
RFstd.to_csv("results/std_ptr-gtex_randomforest_codusnorm_to_UPDOWNproteins.csv")
