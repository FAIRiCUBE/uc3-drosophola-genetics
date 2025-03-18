#!/bin/bash

# Output file
echo "Source, RDA1, RDA2, RDA3" > merged_loadings.csv

# Iterate over directories with the pattern "permut*"
for dir in permut*/; do
    # Check if loadings.csv exists in the directory
    if [[ -f "$dir/loadings.csv" ]]; then
        # Extract the directory name as the source identifier
        source_name=$(basename "$dir")
        # Append the data to the merged file with new headers
        tail -n +2 "$dir/loadings.csv" | awk -v src="$source_name" '{print src "," $0}' >> merged_loadings.csv
    fi
done

echo "Merging complete. Output saved in merged_loadings.csv"




## do the same with pvalues 

# Output file
echo "Source, pvalue" > merged_pvalues.csv

# Iterate over directories with the pattern "permut*"
for dir in permut*/; do
    # Check if loadings.csv exists in the directory
    if [[ -f "$dir/rdadapt_pvalues.csv" ]]; then
        # Extract the directory name as the source identifier
        source_name=$(basename "$dir")
        # Append the data to the merged file with new headers
        tail -n +2 "$dir/rdadapt_pvalues.csv" | awk -v src="$source_name" '{print src "," $0}' >> merged_pvalues.csv
    fi
done

echo "Merging complete. Output saved in merged_loadings.csv"

#### other format

# Initialize output file
echo -n "Index" > all_pvalues.csv

# Collect all directories matching "permut*"
dirs=(permut*/)

# Add directory names as column headers
for dir in "${dirs[@]}"; do
    if [[ -f "$dir/rdadapt_pvalues.csv" ]]; then
        echo -n ",$(basename "$dir")" >> all_pvalues.csv
    fi
done

echo "" >> all_pvalues.csv

# Get the maximum number of rows in any file
max_rows=0
for dir in "${dirs[@]}"; do
    if [[ -f "$dir/rdadapt_pvalues.csv" ]]; then
        rows=$(wc -l < "$dir/rdadapt_pvalues.csv")
        (( rows > max_rows )) && max_rows=$rows
    fi
done

# Merge columns by row index
for ((i=2; i<=max_rows; i++)); do
    echo -n "$((i-1))" >> all_pvalues.csv
    for dir in "${dirs[@]}"; do
        if [[ -f "$dir/rdadapt_pvalues.csv" ]]; then
            value=$(sed -n "${i}p" "$dir/rdadapt_pvalues.csv")
            echo -n ",$value" >> all_pvalues.csv
        else
            echo -n "," >> all_pvalues.csv
        fi
    done
    echo "" >> all_pvalues.csv
done

echo "Merging complete. Output saved in merged_loadings.csv"
