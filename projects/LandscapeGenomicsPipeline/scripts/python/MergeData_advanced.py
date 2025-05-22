import pandas as pd
from optparse import OptionParser, OptionGroup

# Help message and options
usage = "python %prog --metadata file --samplenames file --output filepath [--variable file1 col1 [orig_col1] file2 col1 [orig_col2] file3 col9 [orig_col3]]"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "<put description here>")

# Add options
parser.add_option("--variable1", dest="VAR1", help="Input Variable file and column")
parser.add_option("--col1", dest="COL1", help="Input Variable file and column")
parser.add_option("--variable2", dest="VAR2", help="Input Variable file and column")
parser.add_option("--col2", dest="COL2", help="Input Variable file and column")
parser.add_option("--output", dest="OUT", help="Output file path")
parser.add_option("--samplenames", dest="SAMPLES", help="Samplenames file")
parser.add_option("--metadata", dest="META", help="Metadata file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

# Read the metadata and samples file
dfMeta = pd.read_csv(options.META)
dfSamplelist = pd.read_csv(options.SAMPLES, header=None)
dfSamplelist.columns = ['sampleId']
dfRASDA = pd.read_csv(options.VAR1)
dfRASDA['sampleId'] = dfRASDA['sampleId'].str.strip()
dfWorlclim = pd.read_csv(options.VAR2)
dfRASDA['sampleId'] = dfRASDA['sampleId'].str.strip()

#dfWC=pd.read_csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/dest.worldclim.csv")
# Clean up whitespace in the sample IDs
dfMeta['sampleId'] = dfMeta['sampleId'].str.strip()
dfSamplelist['sampleId'] = dfSamplelist['sampleId'].str.strip()

# Initialize the merged DataFrame with metadata and samples
merged_df = pd.merge(dfMeta, dfSamplelist, on='sampleId', how='inner')
dfWorlclim.rename(columns={'sampleId':'sampleId_orig'}, inplace=True)

merge2 = pd.merge(merged_df, dfWorlclim, on='sampleId_orig', how='inner')
#print(merged_df)
merge3 = pd.merge(merge2, dfRASDA, on='sampleId', how='inner')
# Process each variable file and column

    #var_file, var_col, *orig_col = var_info
    #orig_col = orig_col[0] if orig_col else 'sampleId'
    #print(orig_col)
    #dfVariable = pd.read_csv(var_file)
    #dfVariable[orig_col] = dfVariable[orig_col].str.strip()
    ## Merge based on the specified column
selected_columns = ['sampleId', 'lat', 'long','nFlies','bio1', varname]
selected_columns = ['sampleId', 'lat', 'long','nFlies','bio1', 'near_surface_air_temperature']
result_df = merge3[selected_columns]

# Display or save the result
#result_df = merged_df.drop_duplicates()  # Drop duplicates in case of multiple matches
result_df.to_csv("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1/metadata.csv", index=False)
