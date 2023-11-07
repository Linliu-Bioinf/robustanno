#!/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
import random

"""
__author__ == "liulin4"
__date__ == 2022-05-20
"""

def rotate_theta(adata, theta = 90):
    sin_theta = np.sin(theta * np.pi/180)
    cos_theta = np.cos(theta * np.pi/180)
    rotate_mat = np.array([[cos_theta, - sin_theta], [sin_theta, cos_theta]])
    sp_id = np.dot(rotate_mat, adata.obsm["spatial"].T)
    adata.obsm["spatial"] = sp_id.T
    return adata


def show_cm(cm_df, normalize = True, title = "Confusion matrix", cmap = plt.cm.Reds, fig_size = [8, 8], save = None):
    classes = cm_df.index.tolist()
    cm = cm_df.values
    if normalize:
        cm = cm.astype('float') / cm.sum(axis = 1)[:, np.newaxis]
    else:
        print("Without normalization!")
    fig, ax = plt.subplots(figsize = fig_size)
    im = ax.imshow(cm, interpolation = 'nearest', cmap = cmap)
    ax.figure.colorbar(im, ax = ax)
    ax.set(xticks=np.arange(cm.shape[1]), yticks = np.arange(cm.shape[0]), xticklabels = classes, yticklabels = classes, title = title, ylabel = "True label", xlabel = "Predicted label")
    plt.setp(ax.get_xticklabels(), rotation = 45, ha = "right", rotation_mode = "anchor")
    fmt = '.1f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt), ha = "center", va = "center", \
color = "white" if cm[i, j] > thresh else 'black', fontsize = 6)
    fig.tight_layout()
    plt.xticks(fontsize= 8)
    plt.yticks(fontsize = 8)
    plt.xlim(-0.5, cm.shape[0] -0.5)
    plt.ylim(cm.shape[1]-0.5, -0.5)
    if save:
        plt.savefig(save, dpi = 300)
        print("Figure has been saved as {}".format(save))


def show_coreCells(adata, scores = "celltype_scores", save = None):
    plt.figure(figsize = (4, 3))
    N, bins, patches = plt.hist(adata.obs[scores], 20, rwidth = 1, alpha = 0.5)
    colors = [(0, 0, 0), "b", "green", "red"]
    divisions = [range(5), range(5, 15), range(15, 20)]
    labels = ["Low quality", "Intermediate cells", "Core cells"]
    for d in divisions:
        patches[list(d)[0]].set_label(labels[divisions.index(d)])
        for i in d:
            patches[i].set_color(colors[divisions.index(d)])
    plt.title("Annotation scores")
    plt.xlabel("Scores")
    plt.ylabel("Cell number")
    plt.legend()
    if save:
        plt.savefig(save)
        print("Scoring has been saved as {}".format(save))
    thre1 = bins[5]
    thre2 = bins[15]
    low_num = adata.obs[scores][adata.obs[scores] <= thre1].shape[0]
    core_num = adata.obs[scores][adata.obs[scores] >= thre2].shape[0]
    inter_num = adata.obs[scores].shape[0] - low_num - core_num
    print("Low quality threshold score is {}".format(thre1))
    print("low celss: {}, inter cells: {}, core cells: {}".format(low_num, inter_num, core_num))


def show_individualScore(adata, key = "celltype_pred", scores = "celltype_scores", a_type = None, save = None):
    if a_type:
        ad = adata[adata.obs[key] == a_type]
        sns.histplot(ad.obs[scores])
        if save:
            plt.savefig(save, dpi = 200)
            print("Type {} scoring is saved as {}".format(a_type, save))       
    else:
        if adata.obs[key].unique() >= 8:
            random_n = 8
        else:
            random_n = adata.obs[key].unique()

        types = random.sample(list(adata.obs[key].unique()), random_n)
        a = 4
        b = 2
        fig, axes = plt.subplots(nrows = a, ncols = b, constrained_layout = True, figsize = (8, 4))
        for nn, ax in enumerate(axes.flat):
            i = types[nn]
            ad = adata[adata.obs[key] == i]
            sns.histplot(ad.obs[scores], ax = ax)
            ax.set(xlabel=' ', ylabel=' ')
            ax.set_title("{}".format(i), fontsize = 8)
        fig.supxlabel("Scores")
        fig.supylabel("Cell number")
        fig.suptitle("Annotation scores")
        plt.subplots_adjust(hspace=0.4, wspace = 0.4)
        if save:
            plt.savefig(save, dpi = 200)
            print("Random {} types scoring has been saved as {}".format(random_n, save))



