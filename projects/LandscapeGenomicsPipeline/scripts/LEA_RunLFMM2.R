
args <- commandArgs(TRUE)
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args){
  cat("Example: -working directory ")}


#BiocManager::install(version = '3.16')
#BiocManager::install("LEA")
#install.packages('gtools')
library(LEA)
library(readr)
library(gtools)

dir.create(args[1])
setwd(args[1])
fol <- paste(args[1])

DATA_full <- read.table(args[2], h=T)

#DATA_full <- read.table(args[2], h=T, sep="\t")
SUB <- DATA_full[3:length(DATA_full)]
samples <- colnames(SUB)
#print(colnames(SUB))
#print(SUB[,1])
metadata <- read.csv(args[3])
#print(SUB)
samples <- colnames(SUB)
samples <- gsub("\\.", "-", samples)
metadata <- read.csv(args[3])
#metadata2 <- read.csv("/media/inter/ssteindl/FC/DEST_freeze1/populationInfo/worldClim/dest.worldclim.csv")
s_data <- metadata[metadata$sampleId %in% samples,]

#write lfmm genotypes 
tr <- t(as.matrix(SUB))
tr <- round(tr,4)

### one #
write.lfmm(tr,"genotypes.lfmm")
#
#### NEWLY ADDED ####
#pc = pca("genotypes.lfmm")
#tw = tracy.widom(pc)
#a=stars.pval(tw$pvalues)
#plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)
#nK<-quick.elbow(tw$eigenvalues, low=0.02, max.pc = 0.9)
#nk=7
#genotype=lfmm2geno("genotypes.lfmm")

print(args[4])
#write environmental lfmm with K number of latent factos and r number of repetitions
var=grep(args[4], colnames(s_data))
varname=paste(args[4])
##K=args[5]
r=args[5]
####this factor needs to be possible as parameter
##factor=s_data[,3]
f <- s_data[,var]
#print(factor)
write.env(f, "gradients.env")

print('CALCULATING LFMM2')
#
mod2 <- lfmm2("genotypes.lfmm", 
                "gradients.env", 
                K = 7)
pv <- lfmm2.test(object = mod2, input = "genotypes.lfmm", env = "gradients.env", linear = TRUE)
#plot(-log10(pv$pvalues), col = "grey", cex = .4, pch = 19)
markerpos<- as.data.frame(cbind(DATA_full$Chr, DATA_full$Pos))
markerpos <- cbind(markerpos,pv$pvalues)

write.csv(markerpos, paste(varname,"_LEA_pvals.csv", sep=""))
