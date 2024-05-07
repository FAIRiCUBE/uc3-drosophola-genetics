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
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numpy.random import randint
from numpy.random import seed
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

write_file=True

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/data/raw/"

out_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/documentation/"

# in_file_name="Europe_50kMutations_5perc_missing.csv"
in_file_name="Europe_50kMutations_gap5perc_expdist.csv"
in_class_name="Europe_50kMutations_t_class5.csv"
out_file_name="Europe_50kMutations_gap5perc_expdist_kmeansfill5.csv"

# in_file_name="North_America_50kMutations_5perc_missing.csv"
# in_file_name="North_America_50kMutations_gap5perc_expdist.csv"
# in_class_name="North_America_50kMutations_t_class5.csv"
# out_file_name="North_America_50kMutations_gap5perc_expdist_kmeansfill5.csv"


df = pd.read_csv(in_base_name+in_file_name,sep=',')
df_class = pd.read_csv(in_base_name+in_class_name,sep=',')
df_gf=df.copy()


[dim_x, dim_y]=df.shape

for ii in range(2,dim_y):
    index_name=df.iloc[:,ii].name
    index_class=np.squeeze(df_class[df_class["Unnamed: 0"]==index_name]["Class"].values)
    class_names=df_class[df_class["Class"]==index_class]["Unnamed: 0"]
    numbers_of_nan=np.size(df.loc[df[index_name].isna(),index_name])
    print("Numbers of NaN found in column ",index_name," : ",numbers_of_nan)
    if numbers_of_nan>0:
        df_gf.loc[df_gf[index_name].isna(),index_name]=df.loc[df_gf[index_name].isna(),class_names].mean(axis=1)
    

# # seed random number generator
# seed(1)
# # generate some integers
# number_of_gaps=3400
# values_pop = randint(2, dim_y, number_of_gaps)
# values_locus = randint(0, dim_x, number_of_gaps)
# values_locus_length=randint(1, 500, number_of_gaps)

# for ii in range(np.size(values_pop)):
#     df_gap.iloc[values_locus[ii]:values_locus[ii]+values_locus_length[ii],values_pop[ii]]=np.NaN

print("Number of NaN in input file : ", df.isna().sum().sum())
print("Number of NaN in output file : ", df_gf.isna().sum().sum())
# print("Percentage of NaN in ouput file : ", df_gap.isna().sum().sum()/df_gap.size*100)

if write_file:
    print("Write CSV output file ",in_base_name+out_file_name)
    df_gf.to_csv(in_base_name+out_file_name,index=False)
    
    
##### stop monitoring of compute resources
# it is very important to use program_path = __file__
shape = df.shape
measurer.end(tracker=tracker,
              shape=shape,
              libraries=[v.__name__ for k, v in globals().items() if type(v) is ModuleType and not k.startswith('__')],
              data_path=data_path,
              program_path=__file__,
              variables=locals(),
              csv_file=out_base_name+filename_mon+'.csv')
