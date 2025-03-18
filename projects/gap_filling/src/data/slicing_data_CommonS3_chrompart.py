#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 11:13:26 2024

@author: sjet
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import boto3
from io import BytesIO
from measurer import Measurer
from types import ModuleType
import os

try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

# Read AWS credentials from environment variables
aws_access_key_id = os.getenv('S3_FAIRICUBE_STORAGE_KEY')
aws_secret_access_key = os.getenv('S3_FAIRICUBE_STORAGE_SECRET')
aws_bucket = os.getenv('S3_FAIRICUBE_STORAGE_BUCKET')

# aws_access_key_id = os.getenv('S3_USER_STORAGE_KEY')
# aws_secret_access_key = os.getenv('S3_USER_STORAGE_SECRET')
# aws_bucket = os.getenv('S3_USER_STORAGE_BUCKET')

# Initialize boto3 client with credentials
s3 = boto3.client(
    's3',
    aws_access_key_id=aws_access_key_id,
    aws_secret_access_key=aws_secret_access_key
)

##### start monitoring of compute resources
out_file_name_mon = os.path.basename(sys.argv[0])
out_base_name_mon = os.getcwd()

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

write_file = True

input_prefix = "genomics_data/raw/filtered/"
output_prefix = "genomics_data/windows/chromosome_section/"

in_file_name_raw1 = "Europe_genomewide_freq.csv"
in_file_name_raw2 = "Europe_genomewide_weight.csv"
in_file_name_index = "Europe_chromatin"

out_file_name1 = "Europe_chromopart_freq"
out_file_name2 = "Europe_chromopart_weight"

# Function to read data from S3
def read_data_from_s3(bucket_name, file_key, se=','):
    obj = s3.get_object(Bucket=bucket_name, Key=file_key)
    data = obj['Body'].read()
    df = pd.read_csv(BytesIO(data), sep=se)
    return df

# Function to write data to S3
def write_data_to_s3(bucket_name, file_key, df):
    csv_buffer = BytesIO()
    df.to_csv(csv_buffer, index=False)
    s3.put_object(Bucket=bucket_name, Key=file_key, Body=csv_buffer.getvalue())

# Ensure the output folder exists in S3
def ensure_folder_exists(bucket_name, folder_key):
    if not folder_key.endswith('/'):
        folder_key += '/'
    s3.put_object(Bucket=bucket_name, Key=folder_key, Body='')

# Ensure the output folder exists
ensure_folder_exists(aws_bucket, output_prefix)

print("Loading input files")

df_raw1 = read_data_from_s3(aws_bucket, input_prefix + in_file_name_raw1, "\t")
df_raw2 = read_data_from_s3(aws_bucket, input_prefix + in_file_name_raw2, "\t")
df_ind = read_data_from_s3(aws_bucket, input_prefix + in_file_name_index + ".csv")


ind_uniq = df_ind["Type"].unique()

print(".... done")
#ind_uniq=ind_uniq[1]

for ii in ind_uniq:
    print("Slicing according to chromosome part", ii)

    chrompart_pos=df_ind[df_ind["Type"] == ii].Pos
    chrompart_pos=chrompart_pos.tolist()

    mask=df_raw1["POS"].isin(chrompart_pos)
    # slice through the raw frequency file
    df_slice1 = df_raw1[mask]
    # slice through the raw weight file
    df_slice2 = df_raw2[mask]

    [dim_x, dim_y] = df_slice1.shape
    # sort according to POS number within sliced data
    df_slice1 = df_slice1.sort_values(by=['POS'])
    df_slice2 = df_slice2.sort_values(by=['POS'])

    print("Number of total samples input file : ", df_slice1.shape[0] * df_slice1.shape[1])
    print("Number of NaN in input file : ", df_slice1.isna().sum().sum())

    if write_file:
        print("Write CSV output file ",[aws_bucket, output_prefix + out_file_name1 + "_" + str(ii) + ".csv"]," for index ", ii, " to S3")
        write_data_to_s3(aws_bucket, output_prefix + out_file_name1 + "_" + str(ii) + ".csv", df_slice1)
        write_data_to_s3(aws_bucket, output_prefix + out_file_name2 + "_" + str(ii) + ".csv", df_slice2)

    print("Creating gap distribution")
    # seed random number generator
    np.random.seed(1)
    # generate some integers
    number_of_gaps = int((dim_x * dim_y) * 2.1 / 100)
    values_pop = np.random.randint(2, dim_y, number_of_gaps)
    values_locus = np.random.randint(0, dim_x, number_of_gaps)
    values_locus_length = np.random.exponential(2, number_of_gaps).astype(int) + 1

    print("Applying gap distribution")
    # adding gaps
    for jj in range(np.size(values_pop)):
        df_slice1.iloc[values_locus[jj]:values_locus[jj] + values_locus_length[jj], values_pop[jj]] = np.nan
        df_slice2.iloc[values_locus[jj]:values_locus[jj] + values_locus_length[jj], values_pop[jj]] = np.nan

    print("Number of NaN in output file : ", df_slice1.isna().sum().sum())
    print("Percentage of NaN in output file : ", df_slice1.isna().sum().sum() / df_slice1.size * 100)

    if write_file:
        print("Write CSV output file ",[aws_bucket, output_prefix + out_file_name1 + "_" + str(ii) + "_gap.csv"]," for index ", ii, " to S3")
        #print("Write CSV output file for index ", ii, " to S3 with gaps")
        write_data_to_s3(aws_bucket, output_prefix + out_file_name1 + "_" + str(ii) + "_gap.csv", df_slice1)
        write_data_to_s3(aws_bucket, output_prefix + out_file_name2 + "_" + str(ii) + "_gap.csv", df_slice2)

    del df_slice1, df_slice2

##### stop monitoring of compute resources
# it is very important to use program_path = __file__
print("Write resource monitoring CSV output file to S3 : ",out_base_name_mon+out_file_name_mon,".csv")
measurer.end(tracker=tracker,
              shape=shape,
              libraries=[v.__name__ for k, v in globals().items() if type(v) is ModuleType and not k.startswith('__')],
              data_path=data_path,
              program_path=__file__,
              variables=locals(),
              csv_file=out_base_name_mon+"/"+out_file_name_mon+".csv")
