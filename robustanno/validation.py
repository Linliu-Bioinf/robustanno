#!/bin/env python3

import pandas as pd
import numpy as np
import scanpy as sc
import os 
import matplotlib.pyplot as plt
import anndata
import seaborn as sns
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.neighbors import NearestCentroid
from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.svm import SVC
from sklearn.metrics import ConfusionMatrixDisplay
from collections import defaultdict

"""
__author__ == "liulin4@genomics.cn"
__date__ == "2022-05-20"
"""

def scoring(adata, n_pcs = 20, n_repeats = 100, key = "celltype_pred"):
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata)
    X = adata.obsm["X_pca"][:, :n_pcs]
    if "label" not in adata.obs:
        label, _ = convert_type_to_label(adata, key)
        adata.obs["label"] = np.asarray(label)
    y = np.array(adata.obs["label"].values)
    idx = adata.obs.index
    df = pd.DataFrame(np.zeros((X.shape[0], n_repeats)), columns = ["repeat" + str(i) for i in range(n_repeats)], index = adata.obs.index)
    df.insert(0, "Label", y, allow_duplicates=False)
    df["Label"] = df["Label"].astype(float)
    rkf = RepeatedKFold(n_splits = 5, n_repeats = n_repeats, random_state=2652124)
    for k, (train_idx, test_idx) in enumerate(rkf.split(X, y)):
        clf = NearestCentroid()
        clf.fit(X[train_idx], y[train_idx])
        y_pred = clf.predict(X[test_idx])
        i = int(k/5)
        new_col = pd.DataFrame(y_pred, columns = ["repeat" + str(i)], index = idx[test_idx])
        df.update(new_col)
        if (k+1) % n_repeats == 0:
            print("{} repeat is trained!".format(k+1))
    scores = list()
    for i in range(df.shape[0]):
        scores.append((df.iloc[i][df.iloc[i] == df.iloc[i]["Label"]].shape[0] -1)/n_repeats)
    adata.obs["celltype_scores"] = scores
    return adata


def svm(adata, by = "pca_comp", deg = None, n_pcs = 50, key = "celltype_pred", method = "SVM"):
    label, type_label_dict = convert_type_to_label(adata, key)
    adata.obs["label"] = np.asarray(label)
    data = pd.DataFrame(adata.obs["label"]).groupby("label").sample(n = None, frac = .8, random_state = 1)
    label = np.array(data["label"])
    if by == "pca_comp":
        if "X_pca" not in adata.obsm:
            sc.tl.pca(adata)
        X_train, X_test, y_train, y_test = train_test_split(adata[data.index, :].obsm["X_pca"][:, :n_pcs].toarray(), label, test_size = .2, random_state =42)
    elif by == "deg":
        adata = adata[:, deg]
        X_train, X_test, y_train, y_test = train_test_split(adata[data.index, :].X.toarray(),label, test_size = .2, random_state =42)
    clf = SVC(gamma=0.001)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    recall = metrics.recall_score(y_test, y_pred, average = None)
    precision = metrics.precision_score(y_test, y_pred, average = None, zero_division = "warn")
    cm = metrics.confusion_matrix(y_test, y_pred)
    precision_df = pd.DataFrame({"precision": precision, "recall": recall}, index = np.unique(y_test))
    precision_df.index = [list(type_label_dict.keys())[list(type_label_dict.values()).index(i)] for i in precision_df.index.tolist()]
    df_cm = pd.DataFrame(cm, index = precision_df.index, columns = precision_df.index)
    return precision_df, df_cm


def convert_type_to_label(adata, key = "celltype_pred"):
# The labelled cell type has to be in the adata.obs
    type_label_dict = defaultdict()
    num = 0
    labels = []
    for i, j in enumerate(adata.obs[key]):
        if j not in type_label_dict.keys():
            type_label_dict[j] = num
            labels.append(num)
            num += 1
        else:
            labels.append(type_label_dict[j])
    labels = np.asarray(labels)
    return labels, type_label_dict







