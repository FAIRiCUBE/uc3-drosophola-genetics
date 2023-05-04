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

write_file=True

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/data/raw/"

out_base_name_pic="/home/sjet/repos/uc3-drosophola-genetics/Documentation/"

# in_file_name="Europe_50kMutations_0.05missing.tsv"
in_file_name="Europe_50kMutations.tsv"
out_file_name="Europe_50kMutations_t_class.csv"
# in_file_name="Europe_50kMutations_5perc_missing.tsv"



# def load_data(x):
#     ''' import data either from a gzipped or or uncrompessed file or from STDIN'''

#     if x == "-":
#         y = sys.stdin
#     elif x.endswith(".gz"):
#         y = gzip.open(x, "rt", encoding="latin-1")
#     else:
#         y = open(x, "r", encoding="latin-1")
#     return y

# def check_conversion(x, d_type):
#     xb = np.array(x, dtype= d_type)
#     # print x, xb
#     if np.array_equal(x, xb):
#         return d_type, True
#     return d_type, False

# popslist = []
# for l in load_data(in_base_name+in_file_name):
    
#     if l.startswith("##"):
#         # OUT2.write(l)
#         continue
#     if l.startswith("#"):
#         # OUT2.write(l)
#         a = l.rstrip().split()
#         header = a[9:]
#         continue
    
#     population = l.rstrip().split()[9:]
    
#     # pops = a[9:]
#     # for i in range(len(pops)):
#     popslist.append(population)
#     # # OUT2.write("\t".join(a[:9])+"\t"+"\t".join(popslist)+"\n")
# # OUT1.close()
# # OUT2.close()

# pops=np.array(popslist[:1000])
# pops=np.array(popslist[:])

df = pd.read_csv(in_base_name+in_file_name,sep='\t')

# df.insert(2, "STD", df.iloc[:,2:].std(axis=1), True)

# data_type_a = "float32"
# print(check_conversion(pops[0,:], data_type_a))

#copy of population array witout incomplete records
# pops_float=np.zeros(np.shape(pops))
# pops_sparse=np.zeros(np.shape(pops))
#std of copy of population array witout incomplete records


[dim_x, dim_y]=df.shape


# for ii in range(np.shape(pops)[0]):
#     try:
#        pops_float[ii,:]=pops[ii,:].astype(float)
#        non_zero_index=np.nonzero(pops_float[ii,:])
#        # if ii==240:
#        #      print("dd")
#        if np.shape(non_zero_index)[1]>0:
#            pops_sparse_std[ii]=np.std(pops_float[ii,non_zero_index])
#        else:
#            pops_sparse_std[ii]=0              
            
           
#        # pops_sparse[ii,non_zero_index]=pops_sparse[ii,non_zero_index]-np.mean(pops_sparse[ii,non_zero_index])
#        # if ii==240:
#        #     print("dd")
                              
                                  
#     except:
#        pops_float[ii,:]=np.nan
#        print("conversion error in line ",ii)
# #     np.where(pops[:,ii])
# #     pops.sparse[pops_sparse-]

pops_sparse_std=df.iloc[:,2:].std(axis=1)
print("Number of input genom records ",dim_x)
print("Number of complete genom records ",dim_x-df.isna().sum(axis = 1).sum())


####################################################
# STD of input data
fig, axs = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [3, 1]})

axs[0].plot(pops_sparse_std,'-b')
# df['STD'].plot(ax=axs[0])
# axs[0].plot(pops_sparse_std,'r*')
# plt.yticks(tick_marks2, class_names, rotation=0 )
axs[0].set_xlabel('Allele postition')
axs[0].set_ylabel('STD of Allele Frequency')
# plt.title('Confusion Matrix for Random Forest Model')
plt.show()

axs[1].hist(pops_sparse_std, bins=20)
axs[1].set_xlabel('STD of Allele Frequency')
axs[1].set_ylabel('Frequency')
plt.show()
plt.savefig(out_base_name_pic+'STD_of_Allele_Frequency_hist.png')

#set all std values below the mean threshold to 0
# significant_threshold=np.mean(pops_sparse_std)
significant_threshold=np.quantile(pops_sparse_std,0.995)

pops_sparse_std=np.where(pops_sparse_std<significant_threshold,0, pops_sparse_std)
pops_sparse_index=np.asarray(np.nonzero(pops_sparse_std))

print("Number of significant input genom records",np.count_nonzero(pops_sparse_std))

####################################################
# STD of input data
fig, axs = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [3, 1]})

# plt.plot(pops_sparse_std,'-b')
axs[0].plot(pops_sparse_std,'r*')
# plt.yticks(tick_marks2, class_names, rotation=0 )
axs[0].set_xlabel('Allele postition')
axs[0].set_ylabel('STD of Allele Frequency')
# plt.title('Confusion Matrix for Random Forest Model')
plt.show()

axs[1].hist(pops_sparse_std[np.where(pops_sparse_std>0)], bins=20)
axs[1].set_xlabel('STD of Allele Frequency')
axs[1].set_ylabel('Frequency')
plt.show()
plt.savefig(out_base_name_pic+'STD_of_Allele_Frequency_hist_filt.png')

####################################################
# Min/Max/mean of input data across populations

fig, axs = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [3, 1]})

# plt.plot(pops_sparse_std,'-b')
axs[0].plot(range(0,dim_y-2),df.iloc[:,2:].mean(axis=0),'-b')
# plt.yticks(tick_marks2, class_names, rotation=0 )
axs[0].set_xlabel('Population number')
axs[0].set_ylabel('Mean of Allele Frequency')
# plt.title('Confusion Matrix for Random Forest Model')
plt.show()

axs[1].hist(df.iloc[:,2:].mean(axis=0), bins=20)
axs[1].set_xlabel('Mean of Allele Frequency')
axs[1].set_ylabel('Frequency')
plt.show()
plt.savefig(out_base_name_pic+'Mean_of_Allele_Frequency_hist_pops.png')


####################################################
# Min/Max/mean of input data across allele positions

fig, axs = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [3, 1]})

# plt.plot(pops_sparse_std,'-b')
axs[0].plot(df.iloc[:,2:].mean(axis=1),'g.')
# plt.yticks(tick_marks2, class_names, rotation=0 )
axs[0].set_xlabel('Allele position')
axs[0].set_ylabel('Mean Allele Frequency')
# plt.title('Confusion Matrix for Random Forest Model')
plt.show()

axs[1].hist(df.iloc[:,2:].mean(axis=1), bins=20)
axs[1].set_xlabel('Mean Allele postition')
axs[1].set_ylabel('Frequency')
plt.show()
plt.savefig(out_base_name_pic+'Mean_of_Allele_Frequency_hist.png')

fig, axs = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [3, 1]})

# plt.plot(pops_sparse_std,'-b')
axs[0].plot(df.iloc[:,2:].min(axis=1),'g.')
# plt.yticks(tick_marks2, class_names, rotation=0 )
axs[0].set_xlabel('Allele position')
axs[0].set_ylabel('Minimum Allele Frequency')
# plt.title('Confusion Matrix for Random Forest Model')
plt.show()

axs[1].hist(df.iloc[:,2:].min(axis=1), bins=20)
axs[1].set_xlabel('Minimum Allele postition')
axs[1].set_ylabel('Frequency')
plt.show()
plt.savefig(out_base_name_pic+'Min_of_Allele_Frequency_hist.png')

fig, axs = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [3, 1]})

# plt.plot(pops_sparse_std,'-b')
axs[0].plot(df.iloc[:,2:].max(axis=1),'g.')
# plt.yticks(tick_marks2, class_names, rotation=0 )
axs[0].set_xlabel('Allele position')
axs[0].set_ylabel('Maximum Allele Frequency')
# plt.title('Confusion Matrix for Random Forest Model')
plt.show()

axs[1].hist(df.iloc[:,2:].max(axis=1), bins=20)
axs[1].set_xlabel('Maximum Allele postition')
axs[1].set_ylabel('Frequency')
plt.show()
plt.savefig(out_base_name_pic+'Max_of_Allele_Frequency_hist.png')

# df=pd.DataFrame({"ident" : header})
# data=np.squeeze(np.transpose(pops_float[pops_sparse_index,:]))
# columns=np.squeeze(pops_sparse_index.astype(str))
# df_new=pd.concat([df, pd.DataFrame(data=data, columns=columns)], axis=1)




X = df.copy()
y = df.iloc[:,2:].columns
X.drop(['#CHROM', 'POS'], axis=1, inplace=True)
X_t=X.transpose()


# print("start kmeans with elbow")
max=5 
kmeans = KMeans(n_clusters = max, init = 'k-means++', max_iter = 300, n_init = 10, random_state = 0)
kmeans.fit(X_t)

X_t.insert(0, "Class", kmeans.labels_, True)

# cs = []
# max = 10
# for i in range(1, max):
#     kmeans = KMeans(n_clusters = i, init = 'k-means++', max_iter = 300, n_init = 10, random_state = 0)
#     kmeans.fit(X)
#     cs.append(kmeans.inertia_)

# plt.figure(figsize=(16,7))
# plt.plot(range(1, max), cs)
# plt.title('The Elbow Method')
# plt.xlabel('Number of clusters')
# plt.ylabel('CS')
# plt.show()
# plt.savefig(out_base_name_pic+'kmeasn_elbow_method.png')

# print("start kmeans with silhouette")
# sil = []
# max = 10

# # dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
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
# plt.savefig(out_base_name_pic+'kmeasn_sil_method.png')

if write_file:
    X_t['Class'].to_csv(in_base_name+out_file_name, index=True)
