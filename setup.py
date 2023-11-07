#!/bin/env python3

"""
__date__ = '2022-05-22'
__author__ = 'liulin4@genomics.cn'
"""

import sys
from setuptools import setup, find_packages
from pathlib import Path

_version_ = "1.0.0"


setup(
    name = "robustanno",
    version = _version_,
    author="liulin4",
    author_email = "liulin4@genomics.cn",
    description = "Functions for mouse brain cell type annotation",
    keywords=["Spatial_transcriptomics", "stereo-seq", "bioinformatics"],
    packages = find_packages(),
    classifiers = [
        'Programming Language :: Python :: 3.6'
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
        'Environment :: Console'
    ],
    install_requires = [
    ]
)


