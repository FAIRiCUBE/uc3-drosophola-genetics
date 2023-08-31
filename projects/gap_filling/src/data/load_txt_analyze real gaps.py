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

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/data/external/"
# out_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/documentation/"

# in_file_name_real       ="Europe_PoolSNP.001.50.8Jun2023.norep.gapprofile_byrows.txt"
in_file_name_real       ="North_America_PoolSNP.001.50.8Jun2023.norep.gapprofile_byrows.txt"

df_real = pd.read_csv(in_base_name+in_file_name_real,sep='\t')


# print("Number of NaN in GAP input file : ", df_gap.isna().sum().sum())
# print("Number of NaN in GAPFILL input file : ", df_gapfill.isna().sum().sum())
# print("Percentage of NaN in GAP file : ", df_gap.isna().sum().sum()/df_gap.size*100)
# print("Percentage of NaN in GAPFILL file : ", df_gapfill.isna().sum().sum()/df_gapfill.size*100)

# df_gap_index=np.asarray(df_gap.isnull()).nonzero()
# df_gap=df_gap.fillna(0)

# df_gapfill_index=np.asarray(df_gapfill.isnull()).nonzero()
# df_gapfill=df_gapfill.fillna(0)

# print("Sum of error in GAP file : ", (df_true.iloc[:,2:]-df_gap.iloc[:,2:]).abs().sum().sum())
# print("Sum of error in GAPFILL file : ", (df_true.iloc[:,2:]-df_gapfill.iloc[:,2:]).abs().sum().sum())

# MSE_gap = np.square(np.array((df_true.iloc[:,2:])-np.array(df_gap.iloc[:,2:]))).mean()   
# MSE_gapfill = np.square((np.array(df_true.iloc[:,2:])-np.array(df_gapfill.iloc[:,2:]))).mean()   
# rsme_gap = np.sqrt(MSE_gap)
# rsme_gapfill = np.sqrt(MSE_gapfill)

# # MSE_gaponly = np.square((df_true.iloc[df_gap_index[0],df_gap_index[1]]-df_gap.iloc[df_gap_index[0], df_gap_index[1]])).mean()   

# print("Root Mean Square Error in GAP file : ", rsme_gap.mean())  
# print("Root Mean Square Error in GAPFILL file : ", rsme_gapfill.mean())  
# print(f"Reduction in RMS error by {(rsme_gap.mean()-rsme_gapfill.mean())/rsme_gap.mean()*100:.2f} %")
# # if write_file:
# #     df.to_csv(in_base_name+out_file_name)


# ###########################################
# #############plot all grids, panel of all
# print("#### Plotting file")

# fig, axs = plt.subplots(2, 2, figsize=(15, 12))

fig, axs = plt.subplots(1,2,figsize=(15, 12))
plt.rcParams['axes.grid'] = False

# bin_levels=np.linspace(0.01,100,50)


# im1=axs.plot(np.array(df_real.iloc[0::,2::].mean(axis=1)), color="green")
gaps_perc_perpop=(((np.array(df_real.iloc[0::,1::]).T)*np.array(df_real.iloc[0::,0])).T).sum(axis=0)*100
im1=axs[0].plot(gaps_perc_perpop, color="green")

pop_index_okgap=np.where(gaps_perc_perpop<10)
x_axis=np.linspace(0,np.shape(gaps_perc_perpop)[0]-1,np.shape(gaps_perc_perpop)[0])
im1=axs[0].plot(x_axis[pop_index_okgap],gaps_perc_perpop[pop_index_okgap], color="blue")
# gap_pop_ok=np.mean(np.array(df_real.iloc[0::,1::])[:,pop_index_okgap], axis=2)
gap_pop_ok=np.squeeze(np.array(df_real.iloc[0::,1::])[:,pop_index_okgap])*100

im2=axs[1].plot(gap_pop_ok)


#

# axs[0].title.set_text(in_file_name_gapfill)
axs[1].set_xlim(0, 15)
axs[1].set_ylim(0, 3)
# axs[0].legend(["True Data","GAP Data", "GAP filled Data"])
# axs[1].legend(["True-GAP Data", "True-GAPfilled Data"])

axs[0].set_xlabel('Population number')
axs[0].set_ylabel('Amount of gaps in %')
axs[1].set_xlabel('Gap size')
axs[1].set_ylabel('Count in %')
# axs[0].set_xlabel('Allele frequency')

# # axs[1,0].set_title("True-GAP Dat   
# # axs[0,0].set_title("True Data")
# # axs[0,1].set_title("GAP Data")
# # axs[1,0].set_title("True-GAP Data")
# # axs[1,1].set_title("GAP filled Data")
# # axs[1,0].set_title(city_string_in5)
# # axs[1,1].set_title(city_string_in6)
# # axs[1,2].set_title(city_string_in7)
# # axs[1,3].set_title(city_string_in8)



# plt.show()
# plt.savefig(out_base_name+in_file_name_gapfill+".png")


# print("#### Plotting file done \n")