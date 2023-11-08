#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:18:23 2022

@author: sjet
"""
try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

import sys
import os
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import random
import gzip
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import seaborn as sns
from measurer import Measurer
from types import ModuleType

##### start monitoring of compute resources
data_path = '/'
measurer = Measurer()
tracker = measurer.start(data_path=data_path)
# example -> shape = [5490, 2170]
shape = []
filename_mon=os.path.basename(sys.argv[0])

plt.close('all')
SMALL_SIZE = 10
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

write_file=False

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/data/raw/"
out_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/documentation/"


# in_file_name="Europe_50kMutations_0.05missing.tsv"
# in_file_name="Europe_50kMutations.tsv"
# out_file_name="Europe_50kMutations_t_class20.csv"

in_file_name="North_America_50kMutations.tsv"
out_file_name="North_America_50kMutations_t_class5.csv"



df = pd.read_csv(in_base_name+in_file_name,sep='\t')
[dim_x, dim_y]=df.shape

pops_sparse_std=df.iloc[:,2:].std(axis=1)
print("Number of input genom records ",dim_x)
print("Number of complete genom records ",dim_x-df.isna().sum(axis = 1).sum())




X = df.copy()
y = df.iloc[:,2:].columns
X.drop(['#CHROM', 'POS'], axis=1, inplace=True)
X_t=X.transpose()


# print("start kmeans with elbow")
max=5 
kmeans = KMeans(n_clusters = max, init = 'k-means++', max_iter = 300, n_init = 10, random_state = 0)
kmeans.fit(X_t)

X_t.insert(0, "Class", kmeans.labels_, True)

plt.figure(figsize=(16,7))
plt.hist(kmeans.labels_, bins=5, rwidth=0.8)
plt.xlabel('cluster number')
plt.ylabel('Frequency')
plt.show()
plt.savefig(out_base_name+'NA_kmeans_class_dist.png')

print("start kmeans with elbow")
cs = []
max = 10
for i in range(1, max):
    kmeans = KMeans(n_clusters = i, init = 'k-means++', max_iter = 300, n_init = 10, random_state = 0)
    kmeans.fit(X)
    cs.append(kmeans.inertia_)

plt.figure(figsize=(16,7))
plt.plot(range(1, max), cs)
plt.title('The Elbow Method')
plt.xlabel('Number of clusters')
plt.ylabel('CS')
plt.show()
plt.savefig(out_base_name+'NA_kmeans_elbow_method.png')

# print("start kmeans with silhouette")
# sil = []
# max = 10

# dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
# for k in range(2, max+1):
#   kmeans = KMeans(n_clusters = k).fit(X)
#   labels = kmeans.labels_
#   sil.append(silhouette_score(X, labels, metric = 'euclidean'))

# plt.figure(figsize=(16,7))
# plt.plot(range(1, max), sil)
# plt.title('The Silhouette Method')
# plt.xlabel('Number of clusters')
# plt.ylabel('CS')
# plt.show()
# plt.savefig(out_base_name_pic+'NA_kmeans_sil_method.png')

if write_file:
    print("Write Class CSV output file ",in_base_name+out_file_name)
    X_t['Class'].to_csv(in_base_name+out_file_name, index=True)

##### stop monitoring of compute resources
# it is very important to use program_path = __file__
measurer.end(tracker=tracker,
              shape=shape,
              libraries=[v.__name__ for k, v in globals().items() if type(v) is ModuleType and not k.startswith('__')],
              data_path=data_path,
              program_path=__file__,
              csv_file=out_base_name+filename_mon+'.csv')