import pandas as pd
import argparse

#########################################################   HELP   #########################################################################
description = "Put your description here"
parser = argparse.ArgumentParser(description=description)

#########################################################   CODE   #########################################################################

parser.add_argument("--output", dest="OUT", help="outfile")
parser.add_argument("--samplenames", dest="SAMPLES", help="samplenames file")
parser.add_argument("--metadata", dest="META", help="metadata file")
parser.add_argument("--worldclim", dest="WORLDCL", help="biovariable from metafile")
parser.add_argument("--biovariable", dest="BIO", nargs='+', help="biovariables from metafile")
parser.add_argument("--climate_extra", dest="CLIMATE", nargs='+', help="biovariables from metafile")
parser.add_argument("--climate_var", dest="CLIMATEVAR", nargs='+', help="biovariables from metafile")

options = parser.parse_args()

# Read the CSV files with headers
dfMeta = pd.read_csv(options.META)
dfSamplelist = pd.read_csv(options.SAMPLES, header=None, names=['sampleId'])
dfWorldclim = pd.read_csv(options.WORLDCL)

# Rename SampleId_orig to sampleId in dfWorldclim for merging


# Optional: Read and concatenate CLIMATE files if provided
if options.CLIMATE:
    climate_dfs = []
    for climate_file in options.CLIMATE:
        climate_dfs.append(pd.read_csv(climate_file))
    dfClimate = pd.concat(climate_dfs, ignore_index=True)
else:
    dfClimate = pd.DataFrame()

# Clean up whitespace in the sample IDs
dfMeta['sampleId'] = dfMeta['sampleId'].str.strip()
dfSamplelist['sampleId'] = dfSamplelist['sampleId'].str.strip()
dfWorldclim['sampleId'] = dfWorldclim['sampleId'].str.strip()

biovar = options.BIO

print(dfMeta.columns)
print(dfSamplelist.columns)
print("WORLDCLIM")
print(dfWorldclim.columns)
print(dfWorldclim['sampleId'])

# Merge based on the 'sampleId' column
merged_df = pd.merge(dfMeta, dfSamplelist, on='sampleId', how='inner')
merge2 = pd.merge(merged_df, dfWorldclim, on='sampleId', how='inner')

# Optional: Merge with CLIMATE data if provided
if not dfClimate.empty:
    merge3 = pd.merge(merge2, dfClimate, on='sampleId', how='inner')
else:
    merge3 = merge2

# Selecting specific columns
selected_columns = ['sampleId', 'lat', 'long', 'nFlies']
selected_columns.extend(options.BIO)

# Optional: Extend with CLIMATEVAR columns if provided
if options.CLIMATEVAR:
    selected_columns.extend(options.CLIMATEVAR)

print(selected_columns)

# Ensure selected columns are in the DataFrame
result_df = merge3[selected_columns]
result_df.dropna(subset=selected_columns, inplace=True)

# To save the result to a new CSV file
result_df.to_csv(options.OUT, index=False)
