import pandas as pd
import argparse

#########################################################   HELP   #########################################################################
description = "Put your description here"
parser = argparse.ArgumentParser(description=description)

#########################################################   CODE   #########################################################################

parser.add_argument("--variable", dest="VAR", help="Input Variable file")
parser.add_argument("--output", dest="OUT", help="outfile")
parser.add_argument("--samplenames", dest="SAMPLES", help="samplenames file")
parser.add_argument("--metadata", dest="META", help="metadata file")
parser.add_argument("--worldclim", dest="WORLDCL", help="biovariable from metafile")
parser.add_argument("--biovariable", dest="BIO", help="biovariable from metafile")

options = parser.parse_args()

# Read the first CSV file with header
dfMeta= pd.read_csv(options.META)
#dfMeta=pd.read_csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/dest_v2.samps_8Jun2023.csv")
dfSamplelist = pd.read_csv(options.SAMPLES)
#dfSamplelist = pd.read_csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/northamerica.csv")
dfSamplelist.columns = ['sampleId']
dfVariable = pd.read_csv(options.VAR)
#dfVariable = pd.read_csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/testNSAT/data/NSAT.csv")
dfWorldclim = pd.read_csv(options.WORLDCL)
#dfWorldclim = pd.read_csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/Europe50k/data/dest.worldclim.csv")

# Clean up whitespace in the sample IDs
dfMeta['sampleId'] = dfMeta['sampleId'].str.strip()
dfSamplelist['sampleId'] = dfSamplelist['sampleId'].str.strip()
dfVariable['sampleId'] = dfVariable['sampleId'].str.strip()

biovar=options.BIO
varname=list(dfVariable.columns)[1]
#biovar=dfMeta[options.BIO]

# Merge based on the 'sample Id' column
merged_df = pd.merge(dfMeta, dfSamplelist, on='sampleId', how='inner')
dfWorldclim.rename(columns={'sampleId':'sampleId_orig'}, inplace=True)


merge2 = pd.merge(merged_df, dfWorldclim, on='sampleId_orig', how='inner')
merge3 = pd.merge(merge2, dfVariable, on='sampleId', how='inner')
# Display or save the result

# Selecting specific columns
selected_columns = ['sampleId', 'lat', 'long','nFlies', varname, biovar]
result_df = merge3[selected_columns]

# To save the result to a new CSV file
result_df.to_csv(options.OUT, index=False)
#result_df.to_csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/northamerica_meta.csv", index=False)