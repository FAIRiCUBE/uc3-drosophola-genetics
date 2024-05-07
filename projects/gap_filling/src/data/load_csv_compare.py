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

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/data/raw/"
out_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/documentation/"


in_file_name_true       ="Europe_50kMutations.tsv"
# in_file_name_gap        ="Europe_50kMutations_5perc_missing.csv"
in_file_name_gap        ="Europe_50kMutations_gap5perc_expdist.csv"
# in_file_name_gapfill    ="Europe_50kMutations_gap5perc_expdist_kmeansfill5.csv"
# in_file_name_gapfill    ="Europe_50kMutations_gap5perc_expdist_vaefilled.csv"
in_file_name_gapfill    ="Europe_50kMutations_gap5perc_expdist_empbase.csv"
# in_file_name_gapfill    ="Europe_50kMutations_5perc_missing_empiricalbase.csv"
# in_file_name_gapfill    ="Europe_50kMutations_5perc_missing_maefilled.csv"


# in_file_name_true       ="North_America_50kMutations.tsv"
# in_file_name_gap        ="North_America_50kMutations_gap5perc_expdist.csv"
# # in_file_name_gapfill    ="North_America_50kMutations_gap5perc_expdist_kmeansfill5.csv"
# # in_file_name_gapfill    ="North_America_50kMutations_gap5perc_expdist_vaefilled.csv"
# in_file_name_gapfill    ="North_America_50kMutations_gap5perc_expdist_empbase.csv"
# in_file_name_gapfill    ="North_America_50kMutations_5perc_missing_kmeansfill5.csv"
# in_file_name_gapfill    ="North_America_50kMutations_5perc_missing_maefilled.csv"
# in_file_name_gapfill    ="North_America_50kMutations_5perc_missing_empiricalbase.csv"
# in_file_name_gapfill    ="North_America_50kMutations_gap5perc_expdist_kmeansfill5.csv"


df_true = pd.read_csv(in_base_name+in_file_name_true,sep='\t')
df_gap = pd.read_csv(in_base_name+in_file_name_gap,sep=',')
# df_gapfill = pd.read_csv(in_base_name+in_file_name_gapfill,sep=',')
df_gapfill = pd.read_csv(in_base_name+in_file_name_gapfill,sep='\t')

df_true=df_true.iloc[1::,:]
df_gap=df_gap.iloc[1::,:]
# df_gapfill=df_gapfill.drop(columns="Unnamed: 0")

print("Number of NaN in GAP input file : ", df_gap.isna().sum().sum())
print("Number of NaN in GAPFILL input file : ", df_gapfill.isna().sum().sum())
print("Percentage of NaN in GAP file : ", df_gap.isna().sum().sum()/df_gap.size*100)
print("Percentage of NaN in GAPFILL file : ", df_gapfill.isna().sum().sum()/df_gapfill.size*100)

df_gap_index=np.asarray(df_gap.isnull()).nonzero()
df_gap=df_gap.fillna(0)

df_gapfill_index=np.asarray(df_gapfill.isnull()).nonzero()
df_gapfill=df_gapfill.fillna(0)

print("Sum of error in GAP file : ", (df_true.iloc[:,2:]-df_gap.iloc[:,2:]).abs().sum().sum())
print("Sum of error in GAPFILL file : ", (df_true.iloc[:,2:]-df_gapfill.iloc[:,2:]).abs().sum().sum())

MSE_gap = np.square(np.array((df_true.iloc[:,2:])-np.array(df_gap.iloc[:,2:]))).mean()   
MSE_gapfill = np.square((np.array(df_true.iloc[:,2:])-np.array(df_gapfill.iloc[:,2:]))).mean()   
rsme_gap = np.sqrt(MSE_gap)
rsme_gapfill = np.sqrt(MSE_gapfill)

# MSE_gaponly = np.square((df_true.iloc[df_gap_index[0],df_gap_index[1]]-df_gap.iloc[df_gap_index[0], df_gap_index[1]])).mean()   

print("Root Mean Square Error in GAP file : ", rsme_gap.mean())  
print("Root Mean Square Error in GAPFILL file : ", rsme_gapfill.mean())  
print(f"Reduction in RMS error by {(rsme_gap.mean()-rsme_gapfill.mean())/rsme_gap.mean()*100:.2f} %")
# if write_file:
#     df.to_csv(in_base_name+out_file_name)


###########################################
#############plot all grids, panel of all
print("#### Plotting file")

# fig, axs = plt.subplots(2, 2, figsize=(15, 12))

fig, axs = plt.subplots(2,1,figsize=(15, 12))
# plt.rcParams['axes.grid'] = False

bin_levels=np.linspace(0.01,1,50)

# im1=axs[0,0].hist(np.array(df_true.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step")
# im1=axs[0,1].hist(np.array(df_gap.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step")
# im1=axs[1,0].hist(np.array(df_true.iloc[0::,2::]).flatten(order='C')-
#                   np.array(df_gap.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step")
# im1=axs[1,1].hist(np.array(df_gapfill.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step")

im1=axs[0].hist(np.array(df_true.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step", color="green")
im2=axs[0].hist(np.array(df_gap.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step", color="red")
im3=axs[0].hist(np.array(df_gapfill.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step", color="blue")

# twinaxs = axs.twinx()
im3=axs[1].hist(np.array(df_true.iloc[0::,2::]).flatten(order='C')-
                  np.array(df_gap.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step", color="magenta")


im5=axs[1].hist(np.array(df_true.iloc[0::,2::]).flatten(order='C')-
                  np.array(df_gapfill.iloc[0::,2::]).flatten(order='C'), bins=bin_levels,histtype="step", color="cyan")


axs[0].title.set_text(in_file_name_gapfill)
# axs[1].set_ylim(0, 50000)
axs[0].legend(["True Data","GAP Data", "GAP filled Data"])
axs[1].legend(["True-GAP Data", "True-GAPfilled Data"])

axs[0].set_xlabel('Allele frequency')
axs[0].set_ylabel('Count of allele-frequncies')
axs[1].set_xlabel('Allele frequency')
axs[1].set_ylabel('Count of allele-frequencies')
axs[0].set_xlabel('Allele frequency')

# axs[1,0].set_title("True-GAP Dat   
# axs[0,0].set_title("True Data")
# axs[0,1].set_title("GAP Data")
# axs[1,0].set_title("True-GAP Data")
# axs[1,1].set_title("GAP filled Data")
# axs[1,0].set_title(city_string_in5)
# axs[1,1].set_title(city_string_in6)
# axs[1,2].set_title(city_string_in7)
# axs[1,3].set_title(city_string_in8)



plt.show()
plt.savefig(out_base_name+in_file_name_gapfill+".png")


print("#### Plotting file done \n")