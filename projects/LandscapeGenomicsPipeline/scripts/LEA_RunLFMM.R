##This scripts purpose will be applying LEA LFMM to a VCF and parallelize workflow via manual z-score pulling

args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args){
  cat("Example: -working directory ")}


#BiocManager::install(version = '3.16')
BiocManager::install("LEA")
library(LEA)
library(readr)

dir.create(args[1])
setwd(args[1])
fol <- paste(args[1])

#replace with final af.gz that still must be generated after pipeline
DATA_full <- read_table(args[2])

SUB <- DATA_full[3:length(DATA_full)]
samples <- colnames(SUB)
metadata <- read.csv(args[3])
#metadata2 <- read.csv("/media/inter/ssteindl/FC/DEST_freeze1/populationInfo/worldClim/dest.worldclim.csv")
s_data <- metadata[metadata$sampleId %in% samples,]
print(s_data)
#write lfmm genotypes 
tr <- t(as.matrix(SUB))
write.lfmm(tr, "genotypes.lfmm")

#write environmental lfmm with K number of latent factos and r number of repetitions
var=grep(args[4], colnames(s_data))
K=args[5]
r=args[6]
###this factor needs to be possible as parameter
#factor=s_data[,3]
factor=s_data[,var]
write.env(factor, "gradients.env")

project = lfmm( "genotypes.lfmm", 
                "gradients.env", 
                K = 3, 
                repetitions = 1, 
                project = "new")
