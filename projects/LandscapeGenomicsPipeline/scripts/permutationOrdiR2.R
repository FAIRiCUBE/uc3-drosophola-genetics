#### permutation OrdiR2 step

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
library(robust)
library(qvalue)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(parallel)

#outdir="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDAResearchPlan"
outdir=args[1]
envdata=args[2]
permut=as.numeric(args[3])

dir.create(paste0(outdir, "/ordiR2step"))

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

RDA0 <- rda(AllFreq ~ 1,  Env) 
print("PERFORMING RDA_FULL")
###   2.b) RDA FULL
RDAfull <- rda(AllFreq ~ .,  Env)


print("PERFORMING ORDIR2STEP")
ModFILE <- paste0(outdir, "/ordiR2step/ordiR2step",permut,".rds")

if (file.exists(ModFILE)) {
  # Load the file
  mod <- readRDS(file=ModFILE)
  message("File ordiR2step_mod loaded successfully.")
} else {
  # Generate the file
  mod <- ordiR2step(RDA0, RDAfull, Pin = 0.01, R2permutations = 1000, R2scope = T, parallel = 20)
  saveRDS(mod, file=ModFILE)
  message("File generated and saved successfully.")
}

best_variables <- all.vars(formula(mod))
write.csv(best_variables, file=paste0(outdir, "/ordiR2step/best_vars_",permut,".csv") )
