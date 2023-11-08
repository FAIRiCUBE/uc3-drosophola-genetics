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
from numpy.random import exponential
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

    write_switch=False
plot_switch=True

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/data/raw/"

out_base_name_pic="/home/sjet/repos/uc3-drosophola-genetics/Documentation/"

# in_file_name="Europe_50kMutations_0.05missing.tsv"
in_file_name="Europe_50kMutations.tsv"
out_file_name="Europe_50kMutations_gap5perc_expdist.csv"
# in_file_name="North_America_50kMutations.tsv"
# out_file_name="North_America_50kMutations_gap5perc_expdist.csv" 


df = pd.read_csv(in_base_name+in_file_name,sep='\t')
df_gap=df.copy()

[dim_x, dim_y]=df.shape

# seed random number generator
seed(1)
# generate some integers
number_of_gaps=338000 #Europe
# number_of_gaps=125000 #North America
values_pop = randint(2, dim_y, number_of_gaps)
values_locus = randint(0, dim_x, number_of_gaps)
values_locus_length=exponential(2, number_of_gaps).astype(int)+1

for ii in range(np.size(values_pop)):
    df_gap.iloc[values_locus[ii]:values_locus[ii]+values_locus_length[ii],values_pop[ii]]=np.NaN

print("Number of total samples input file : ", df.shape[0]*df.shape[1])
print("Number of NaN in input file : ", df.isna().sum().sum())
print("Number of NaN in output file : ", df_gap.isna().sum().sum())
print("Percentage of NaN in ouput file : ", df_gap.isna().sum().sum()/df_gap.size*100)

if write_switch:
    df_gap.to_csv(in_base_name+out_file_name, index=False)
    
if plot_switch:
    print("#### Plotting file")
    fig = plt.figure(figsize=(12, 6))
    plt.hist(values_locus_length,bins=np.linspace(1, 15,15),histtype=u'step')
