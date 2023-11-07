#!/bin/env python3
import scanpy as sc
import numpy as np
import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import os
from sklearn.model_selection import RepeatedKFold

"""
__date__ == "2022-05-22"
__author__ == "liulin4@genomics.cn"
"""

def iterateClustering(adata, n_neighbors = 15, n_pcs = 30, r = 0.2, n_repeats = 100, save = None):
    coo_mtx, cell_prob = subsampling(adata, n_neighbors = n_neighbors, n_pcs = n_pcs, r = r, n_repeats = n_repeats)
    idx = np.random.choice(cell_prob.shape[0], size = int(cell_prob.shape[0]*0.1))
    show_mtx = cell_prob[idx, :][:, idx]
    plt.rcParams["figure.dpi"] = 300
    plt.rcParams["figure.figsize"] = [5, 5]
    ax = sns.clustermap(cell_prob)
    if save:    
        plt.savefig(save)
        print("Probability matrix has been saved as {}".format(save))
    return cell_prob


def filteringCluster(adata, cell_prob, n_neighbors = 15, n_pcs = 30, r =0.2, vali_thred = 0.3, save = None):
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, n_pcs = n_pcs)
    sc.tl.leiden(adata, resolution = r)

    low_quality = []
    for i in adata.obs["leiden"].cat.categories:
        print(f"Finetuning cluster {i}!")
        cluster_ids = adata.obs.index[adata.obs["leiden"] == i]
        cluster_idx = np.where(adata.obs.index.isin(cluster_ids))
        prob_mtx = cell_prob[cluster_idx[0], :]
        prob_mtx = prob_mtx[:, cluster_idx[0]]
        idx = np.random.choice(prob_mtx.shape[0], size = int(prob_mtx.shape[0]*0.1))
        show_mtx = prob_mtx[idx, :][:, idx]
        vali_idx = np.where(np.median(prob_mtx, axis = 1) >= vali_thred)[0]
        vali_prob = prob_mtx[vali_idx,:][:, vali_idx]
        idx = np.random.choice(vali_prob.shape[0], size = int(vali_prob.shape[0]*0.1))
        show_mtx = vali_prob[idx, :][:, idx]
        low_quality.extend(cluster_ids[np.where(np.median(prob_mtx, axis = 1) <= vali_thred)[0]])
    tmpdf = pd.DataFrame(adata.obs["leiden"][adata.obs.index.isin(low_quality)])
    for i in tmpdf["leiden"].unique():
        tmpdf["leiden"].replace(i, "low_quality", inplace = True)
    df = pd.DataFrame(adata.obs["leiden"])
    df.update(tmpdf)
    adata.obs["leiden_fined"] = df
    return adata


def subsampling(adata, n_neighbors = 15, n_pcs = 30, r = 0.2, n_repeats = 20, n_split = 5):
    y = np.arange(0, adata.obs.shape[0])
    coo_mtx = np.zeros((adata.obs.shape[0], )*2)
    rkf = RepeatedKFold(n_splits = n_split, n_repeats = n_repeats)
    for epoch, (train_idx, test_idx) in enumerate(rkf.split(y)):
        print(f"Resampling {epoch + 1} iterations!")
        sample_ids = adata.obs.index[y[train_idx]]
        ad = adata[sample_ids,]
        sc.tl.pca(ad)
        sc.pp.neighbors(ad, n_neighbors = n_neighbors, n_pcs = n_pcs)
        sc.tl.leiden(ad, key_added = "leiden", resolution = r)
        for i in ad.obs["leiden"].cat.categories:
            cluster_ids = ad.obs.index[ad.obs["leiden"] == i]
            cluster_idx = np.where(adata.obs.index.isin(cluster_ids))
            for m, n in enumerate(cluster_idx[0][0:]):
                coo_mtx[n, n] += 1
                for j, k in enumerate(cluster_idx[0][m+1:]):
                    coo_mtx[n, k] += 1
                    coo_mtx[k, n] += 1
    cell_prob = coo_mtx / (epoch+1)
    return coo_mtx, cell_prob

