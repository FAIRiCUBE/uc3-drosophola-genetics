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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numpy.random import randint
from numpy.random import seed


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

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/data/raw/"


# in_file_name_true       ="Europe_50kMutations.tsv"
# in_file_name_gap        ="Europe_50kMutations_5perc_missing.csv"
# in_file_name_gapfill    ="Europe_50kMutations_5perc_missing_kmeansfill5.csv"
in_file_name_true       ="North_America_50kMutations.tsv"
in_file_name_gap        ="North_America_50kMutations_5perc_missing.csv"
in_file_name_gapfill    ="North_America_50kMutations_5perc_missing_kmeansfill5.csv"



df_true = pd.read_csv(in_base_name+in_file_name_true,sep='\t')
df_gap = pd.read_csv(in_base_name+in_file_name_gap,sep=',')
df_gapfill = pd.read_csv(in_base_name+in_file_name_gapfill,sep=',')

print("Number of NaN in GAP input file : ", df_gap.isna().sum().sum())
print("Number of NaN in GAPFILL input file : ", df_gapfill.isna().sum().sum())
print("Percentage of NaN in GAP file : ", df_gap.isna().sum().sum()/df_gap.size*100)
print("Percentage of NaN in GAPFILL file : ", df_gapfill.isna().sum().sum()/df_gapfill.size*100)

df_gap=df_gap.fillna(0)

print("Sum of error in GAP file : ", (df_true.iloc[:,2:]-df_gap.iloc[:,2:]).abs().sum().sum())
print("Sum of error in GAPFILL file : ", (df_true.iloc[:,2:]-df_gapfill.iloc[:,2:]).abs().sum().sum())

MSE_gap = np.square((df_true.iloc[:,2:]-df_gap.iloc[:,2:])).mean()   
MSE_gapfill = np.square((df_true.iloc[:,2:]-df_gapfill.iloc[:,2:])).mean()   
rsme_gap = np.sqrt(MSE_gap)
rsme_gapfill = np.sqrt(MSE_gapfill)


print("Root Mean Square Error in GAP file : ", rsme_gap.mean())  
print("Root Mean Square Error in GAPFILL file : ", rsme_gapfill.mean())  
print(f"Reduction in RMS error by {(rsme_gapfill.mean()/rsme_gap.mean())*100:.2f} %")
# if write_file:
#     df.to_csv(in_base_name+out_file_name)
