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
library(gtools)


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


#### NEWLY ADDED ####
pc = pca("genotypes.lfmm", scale = TRUE)
tw = tracy.widom(pc)
a=stars.pval(tw$pvalues)
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
nK<-quick.elbow(tw$eigenvalues, low=0.02, max.pc = 0.9)
genotype=lfmm2geno("genotypes.lfmm")

#### NEWLY ADDED OVER ####

#write environmental lfmm with K number of latent factors and r number of repetitions
var=grep(args[4], colnames(s_data))
#K=args[5]
K=nK
r=args[5]
###this factor needs to be possible as parameter
#factor=s_data[,3]
factor=s_data[,var]

write.env(factor, "gradients.env")

project = lfmm( "genotypes.lfmm", 
                "gradients.env", 
                K = 12, 
                repetitions = 1, 
                project = "new")

##### only relevant when running as standalone in R Studio
#zs = z.scores(obj.lfmm, K=12)
#
##Combine z-scores using the median
#zs.median = apply(zs, MARGIN = 1, median)
#lambda = median(zs.median^2)/0.456
#lambda
#adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
#hist(adj.p.values, col = "red")
#adj.p.values = pchisq(zs.median^2/.85, df = 1, lower = FALSE)
#
##histogram of p-values
#hist(adj.p.values, col = "green")
#L = 500
#q = 0.1
#w = which(sort(adj.p.values) < q * (1:L)/L)
#candidates = order(adj.p.values)[w]
