#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:13:26 2024

@author: sjet
"""
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
from measurer import Measurer
from types import ModuleType
import os

##### start monitoring of compute resources
out_file_name_mon=os.path.basename(sys.argv[0])
out_base_name_mon=os.getcwd()

data_path = '/'
measurer = Measurer()
tracker = measurer.start(data_path=data_path)
# example -> shape = [5490, 2170]
shape = []


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

in_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/makeWindows/data/"
out_base_name="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/data/windows/"
# out_base_name_mon="/home/sjet/repos/uc3-drosophola-genetics/projects/gap_filling/data/windows/"

in_file_name_raw1       ="Europe_genomewide_freq.csv"
in_file_name_raw2       ="Europe_genomewide_weight.csv"
in_file_name_index     ="Europe_chromosomes"

out_file_name1       ="Europe_chromosomes_freq"
out_file_name2       ="Europe_chromosomes_weight"

print("Loading input files")

df_raw1 = pd.read_csv(in_base_name+in_file_name_raw1, sep="\t")
df_raw2 = pd.read_csv(in_base_name+in_file_name_raw2, sep="\t")

df_ind = pd.read_csv(in_base_name+in_file_name_index+".csv")

ind_uniq=df_ind["Chrom"].unique()

print(".... done")

for ii in ind_uniq:
# ii=ind_uniq[0]
    print("Slicing according to chromosome ",ii)
    
    #slice through the raw frequency file
    df_slice1=df_raw1[df_raw1["#CHROM"]==ii]
    #slice through the raw weight file
    df_slice2=df_raw2[df_raw2["#CHROM"]==ii]
    
    #merge sliced files and create a tuple of fre-weight per position & population
    # df_slice= pd.concat([df_slice1,df_slice2]).stack().groupby(level=[0,1]).apply(tuple).unstack() #to slow and running out of memory
    # df_slice = pd.DataFrame({x: zip(df_slice1[x], df_slice2[x]) for x in df_slice1.columns})
    
    #first two columns are now containing tuples of duplicates, will be replaced by single value
    # df_slice1["POS"] = df_slice1["POS"].map(lambda x: x[0]) 
    # df_slice2["POS"] = df_slice2["POS"].map(lambda x: x[0]) 
    # df_slice1["#CHROM"] = df_slice1["#CHROM"].map(lambda x: x[0]) 
    # df_slice2["#CHROM"] = df_slice2["#CHROM"].map(lambda x: x[0]) 
    
    [dim_x, dim_y]=df_slice1.shape
    #sort according to POS number within sliced data
    df_slice1=df_slice1.sort_values(by=['POS'])
    df_slice2=df_slice2.sort_values(by=['POS'])
    
    print("Number of total samples input file : ", df_slice1.shape[0]*df_slice1.shape[1])
    print("Number of NaN in input file : ", df_slice1.isna().sum().sum())
    
    if write_file:
        print("Write CSV output file for index ",ii," to file ",out_base_name+out_file_name1+str(ii)+".csv")
        print("Write CSV output file for index ",ii," to file ",out_base_name+out_file_name2+str(ii)+".csv")
        # df_gap.to_csv(out_base_name+out_file_name+"_chrom"+str(ii)+"_gap.csv", index=True)
        df_slice1.to_csv(out_base_name+out_file_name1+"_chrom"+str(ii)+".csv", index=True)
        df_slice2.to_csv(out_base_name+out_file_name2+"_chrom"+str(ii)+".csv", index=True)
        
    #create a copy of sliced data to insert gaps
    # df_gap=df_slice.copy()
    
    print("Creating gap distribution")
    # seed random number generator
    np.random.seed(1)
    # generate some integers
    #2.1% of number of samples in matrix makes about total of 5% gaps per sliced data
    number_of_gaps=np.int((dim_x*dim_y)*2.1/100) 
    # number_of_gaps=125000 #North America
    values_pop = np.random.randint(2, dim_y, number_of_gaps)
    values_locus = np.random.randint(0, dim_x, number_of_gaps)
    values_locus_length=np.random.exponential(2, number_of_gaps).astype(int)+1
     
    print("Applying gap distribution")
    #adding gaps 
    for jj in range(np.size(values_pop)):
        df_slice1.iloc[values_locus[jj]:values_locus[jj]+values_locus_length[jj],values_pop[jj]]=np.NaN
        df_slice2.iloc[values_locus[jj]:values_locus[jj]+values_locus_length[jj],values_pop[jj]]=np.NaN
    
       
    print("Number of NaN in output file : ", df_slice1.isna().sum().sum())
    print("Percentage of NaN in ouput file : ", df_slice1.isna().sum().sum()/df_slice1.size*100)
        
    if write_file:
        print("Write CSV output file for index ",ii," to file ",out_base_name+out_file_name1+str(ii)+".csv")
        print("Write CSV output file for index ",ii," to file ",out_base_name+out_file_name2+str(ii)+".csv")
        df_slice1.to_csv(out_base_name+out_file_name1+"_chrom"+str(ii)+"_gap.csv", index=True)
        df_slice2.to_csv(out_base_name+out_file_name2+"_chrom"+str(ii)+"_gap.csv", index=True)
        # df_slice.to_csv(out_base_name+out_file_name+"_chrom"+str(ii)+".csv", index=True)
        
    del df_slice1, df_slice2

##### stop monitoring of compute resources
# it is very important to use program_path = __file__
print("Write resource monitoring CSV output file ",out_base_name_mon+out_file_name_mon+".csv")
measurer.end(tracker=tracker,
              shape=shape,
              libraries=[v.__name__ for k, v in globals().items() if type(v) is ModuleType and not k.startswith('__')],
              data_path=data_path,
              program_path=__file__,
              variables=locals(),
              csv_file=out_base_name_mon+"/"+out_file_name_mon+'.csv')