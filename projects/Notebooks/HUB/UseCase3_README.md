## Use Case 3: Landscape Genomics Analysis in FAIRiCUBE

This README provides a simple overview of the Use Case Work. 
### Overview

The Use Case workflow is based on *BASH scripting*, *Python* and *R* for analyzing genomic variant (SNP) data for multiple samples and loci. 
This documentation shall show how to get environmental data relevant for Landscape Genomics Analysis within the FAIRiCUBE structure and consequently allows performing the analysis 

### Third party programs used in the pipeline

* [`NanoFilt`](https://github.com/wdecoster/nanofilt) for read filtering
* [`vcftools`](https://vcftools.sourceforge.net/) for handling genomic data in the variant call format
* [`R`](https://www.r-project.org/) for plotting and statistical analyses


### Execution 

#### SetUp 
[**Notebook 1**: SetUp The Environments and Packages](/user/home/sonjastndl/s3/Notebooks/LandscapeGenomicsSetup.ipynb)

#### Genomic Data 
[**Notebook 2**: Download Genomic Data](s3/Notebooks/GetGenomicData.ipynb)

Local genetic data can be uploaded in the form of a VCF on genotypes or an corresponding allele frequency file.

#### Get Environmental Data  

[**Notebook 3**: Get Environmental Data via Sentinel Hub](s3/ClimateData/DataFromSH.ipynb)

[**Notebook 4**: Get Environmental Data: Alternatives](s3/Notebooks/OOWP.ipynb)

#### Landscape Analysis  
[**Notebook 5**: Run the Landscape Genomics Analysis](s3/Notebooks/Run_LG_Pipeline.ipynb)

## Objectives

These Jupyter Notebooks are outcomes of the Use Case 3 for FAIRiCUBE and document an Environmental Association Analysis Workflow based on Environmental Data obtained and processed within the FAIRiCUBE project and a Landscape Genomics Pipeline constructed within the project.
By following the Workflows of the Notebooks, one can perform an examplatory Analysis on *Drosophila melanogaster* genomic data from DrosEU from 2024. 