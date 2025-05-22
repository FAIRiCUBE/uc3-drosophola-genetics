#install.packages('pegas')
args <- commandArgs(TRUE)
library(dplyr)
library(pegas)
library(ggplot2)
library(raster)
library(ggrepel)
library(LEA)
library(rnaturalearth)
library(rnaturalearthdata)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(ggVennDiagram)
library(cowplot)
library(corrplot)
library(robust)
#BiocManager::install("qvalue")
library(qvalue)
library(ggpubr)
library(geodata)
library('corrr')
#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("FactoMineR")
library(FactoMineR)
#install.packages("factoextra")
library(factoextra)
#library(dplyr)
library(tidyverse)
library(UpSetR)

#outdir="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDAResearchPlan"
outdir=args[1]
envdata=args[2]
permut=as.numeric(args[3])

dir.create(paste0(outdir, "/permut", permut))

# READING ENVIRONMENTAL DATA
#envdata="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/ClimateData/Env_imputed.csv"

Env <- read.csv(envdata, header=TRUE)
rownames(Env) <- Env$sample
Dates <- as.Date(Env$Date_num)
custom_origin <- as.Date(min(Env$Date_num))
Env$Date_num <- as.numeric(as.Date(Dates) - custom_origin)
## Dates <- as.Date(Env$Date_num + custom_origin)
UncertaintyParameters <- c("SMcVSMU", "SMpVSMU", "AImPP", "SMpDNF", "SMcDNF", "SMaDNF", "SMaPercSatU", "AI340P", "AI354P", "AImHP", "CH4mrP", "CLgCBHP", "SO2gVCP", "O3gOVCP", "NO2gAMFtcpk", "NO2gAMFtcp", "HCOHgFVCP", "COtcP", "CLgCTPP", "CLgCTHP", "CLgCOTP", "CLgCFP", "CLgCBPP")
Env <- Env[!colnames(Env) %in% UncertaintyParameters]

Env_backup <- Env


af_file="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Subsampled_fullgenome.final_DP15.af"
af <- read.table(af_file, header=TRUE, na.strings = NaN)
ann="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/data/annotationdata/annotations.txt"
annotations <- read.table(ann, na.strings = NaN)
annotations <- unique(annotations)

colnames(annotations) <- c("Chr", "Pos", "Gene")
rownames(annotations) <- paste0(annotations$Chr,".", annotations$Pos)
loci <- paste(af$Chr, af$Pos, sep=".")
af2 <- as.data.frame(t(af %>% dplyr::select(3:ncol(af))))
colnames(af2) <- loci

AllFreq=af2
### here input args 
permutationRowName <- read.table("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/PermutationMatrix2.rds", sep=" ")[permut][1]
newnames <- as.character(permutationRowName[,1])
rownames(AllFreq) <- gsub("\\.", "-", newnames)


freq_mean <- colMeans(AllFreq)
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]

AllFreq <- AllFreq[rownames(AllFreq) %in% rownames(Env),]
Env <- Env[rownames(Env) %in% rownames(AllFreq),]
Dates <- Env$Date_num
Lat <- Env$lat
Long <- Env$long
Coordinates <- as.data.frame(cbind(Lat, Long))
rownames(Coordinates) <- Env$sample
colnames(Coordinates) <- c("Lat", "Long")
Env <- as.data.frame(Env[,5:ncol(Env)])
Env$Date_num <- Dates

print("Working with Allele Frequencies Means >=0.95 or <=0.05 ")
print(nrow(AllFreq))
print(ncol(AllFreq))

AFsin <- asin(sqrt(AllFreq))

print("SCALING AND CENTERING")
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')


print("FIRST PCA")

pr_eur <- prcomp(Env, scale. = TRUE)
Env.PCA <- princomp(Env)
Env.PCA2 <- pr_eur

Env <- as.data.frame(Env)
AllFreq <- AFsin

####  ORDI2STEP: VARIABLE SELECTION
##done before

bv <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/ordiR2step/best_variables.csv")
###best vars and save
#best_variables <- all.vars(formula(mod))
#best_variables <- setdiff(best_variables, "AllFreq")
best_variables <- bv$x[ bv$x != "Date_num"]

#xo <- ordiplot(mod, scaling = 1, type = "text")

EnvComp <- Env[best_variables]

### PREPARING NEUTRAL SNP AND COORDINATES
neutral_SNPS <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Subsampled_fullgenome_NeutralSNPS_80.tsv"

AllFreq_neutral <- read.table(neutral_SNPS, header = F)
rownames(AllFreq_neutral) <- paste(AllFreq_neutral$V1, AllFreq_neutral$V2, sep=".")

Loci <- rownames(AllFreq_neutral)
AF2_neutral <- AllFreq[,colnames(AllFreq) %in% Loci]

AllFreq_neutral=AF2_neutral[gsub("\\.", "-", rownames(AF2_neutral)) %in% rownames(Env),]

### PERFORMING PCA NEUTRAL

FilePCANeutral <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDA_ResearchPlan_FilterMethod1/populationStructure/pcaNeutral.rds"

pca_neutral <- readRDS(file=FilePCANeutral)

#here first three are enogh to descirbe neutral population genetic structure??
PCs_neutral <- scores(pca_neutral, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = rownames(AllFreq_neutral), PCs_neutral)


print("INTERSECTING COORDINATES AND ENVCOMP AND NEUTRAL PCS")

Variables <- data.frame(Coordinates, EnvComp, PopStruct$PC1, PopStruct$PC2, PopStruct$PC3)
#Variables <- Variables[Variables$sample %in% rownames(AllFreq),]

##only components without lat, long
VariablesN <- Variables[,3:ncol(Variables)]

##included
Factors <- colnames(EnvComp)


### 7.1) Make an RDA without Geography to invest deeper the climatic factors

formula_string <- paste("AllFreq ~ ", paste(Factors, collapse = " + "), " + Condition(PopStruct.PC1 +PopStruct.PC2 + PopStruct.PC3)")
formula_env <- as.formula(formula_string)

RDA_env <- rda(formula_env, Variables)

### 9 IDENTIFY CANDIDATE SNPS INVOLVED IN LOCAL ADAPTATION - APPROACH EXCURSE
Variables <- data.frame(Coordinates, EnvComp, PopStruct$PC1, PopStruct$PC2, PopStruct$PC3)

load.rda <- scores(RDA_env, choices=c(1:3), display="species")  # Species scores for the first three constrained axes
#pvalue <- 0.05
#pvalue <- 0.05/length(rdadapt_env$p.values)
pvalue <- 0.05/ncol(AllFreq)
zscore <- qnorm(1 - pvalue / 2)


#outliers <- function(x,z){
#  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
#  x[x < lims[1] | x > lims[2]]               # locus names in these tails
#}

#cand1 <- outliers(load.rda[,1],zscore) 
#cand2 <- outliers(load.rda[,2],zscore) 
#cand3 <- outliers(load.rda[,3],zscore) 
#

loadings1 <- load.rda[,1]
loadings2 <- load.rda[,2]
loadings3 <- load.rda[,3]

write.csv(load.rda, paste0(outdir, "/permut", permut,"/loadings.csv"), row.names = FALSE)

source("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/src/rdadapt.R")
rdadapt_env<-rdadapt(RDA_env, 2)
write.csv(rdadapt_env$p.values, paste0(outdir, "/permut", permut,"/rdadapt_pvalues.csv"), row.names = FALSE)

