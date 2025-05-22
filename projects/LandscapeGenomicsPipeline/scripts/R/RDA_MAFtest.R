#install.packages('pegas')
library(pegas)
library(ggplot2)
library(raster)
library(ggrepel)
#install.packages("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/packages/rgdal_1.6-7.tar.gz", repos=NULL, type='source')
#library(rgdal)

library(LEA)
#install.packages('rnaturalearth','rnaturalearthdata', 'robust')
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
#library(WMDB)
#library(rgeos)
#install.packages("corrr")
library('corrr')
#install.packages("ggcorrplot")
library(ggcorrplot)
#install.packages("FactoMineR")
library(FactoMineR)
#install.packages("factoextra")
library(factoextra)
#library(dplyr)
library(tidyverse)



# Functions
get_meta_data <- function(csvfile) {
  meta <- read.csv(csvfile, header = TRUE)
  meta.sub <- meta %>% dplyr::select(sampleId, continent, country, province, lat, long)
  return(meta.sub)
}


get_worldclim_data <- function(meta.sub) {
  biod <- worldclim_global(var = "bio", 2.5, "data")
  bio <- raster::extract(biod, cbind(meta.sub$long, meta.sub$lat))
  bio.sub <- as.data.frame(bio)
  return(bio.sub)
}

args <- commandArgs(TRUE)
#af_file<-args[1]
#meta2024<-args[2]
#outdir<-args[4]
#region<-args[5]
#errorperc<- as.numeric(args[6])
#dir.create(outdir)

#meta <- get_meta_data(meta2024)

#af <- read.table(af_file, header=TRUE, na.strings = NaN)
#loci <- paste(af$Chr, af$Pos, sep=".")
#af_sel <- af %>% dplyr::select(3:ncol(af))
#af2 <- as.data.frame(t(af_sel))
#colnames(af2) <- loci
#
#AllFreq=af2
#rownames(AllFreq) <- gsub("\\.", "-", rownames(AllFreq))
#common_elements <- intersect(meta$sampleId, row.names(AllFreq))
## Subset the data frames based on the intersection
#matchedMeta  <- meta[meta$sampleId %in% common_elements, ]


errorperc <- as.numeric(args[1])
outdir <- args[2]
# READING ENVIRONMENTAL DATA
envdata="/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/ClimateData/Env.csv"
Env <- read.csv(envdata, header=TRUE)
rownames(Env) <- Env$sample
Dates <- as.Date(Env$Date_num)
custom_origin <- as.Date(min(Env$Date_num))
Env$Date_num <- as.numeric(as.Date(Dates) - custom_origin)
## Dates <- as.Date(Env$Date_num + custom_origin)
#UncertaintyParameters <- c("SMcVSMU", "SMpVSMU", "AImPP", "SMpDNF", "SMcDNF", "SMaDNF", "SMaPercSatU", "AI340P", "AI354P", "AImHP", "CH4mrP", "CLgCBHP", "SO2gVCP", "O3gOVCP", "NO2gAMFtcpk", "NO2gAMFtcp", "HCOHgFVCP", "COtcP", "CLgCTPP", "CLgCTHP", "CLgCOTP", "CLgCFP", "CLgCBPP")
#Env <- Env[!colnames(Env) %in% UncertaintyParameters]

# READING ALLELE FREQUENCY FILE 
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
rownames(AllFreq) <- gsub("\\.", "-", rownames(AllFreq))
freq_mean <- colMeans(AllFreq)

## WORK

errorval_lower <- 1*errorperc/100
errorval_higher <- 1 - errorval_lower

#AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]
AllFreq <- AllFreq[,-which(freq_mean>=errorval_higher | freq_mean<=errorval_lower)]

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

print("Testing different Allele Frequencies")
print(nrow(AllFreq))
print(ncol(AllFreq))

AFsin <- asin(sqrt(AllFreq))

print("SCALING AND CENTERING")
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')


### PREPARING NEUTRAL SNP AND COORDINATES
neutral_SNPS <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/Subsampled_fullgenome_NeutralSNPS_80.tsv"


AllFreq_neutral <- read.table(neutral_SNPS, header = F)
rownames(AllFreq_neutral) <- paste(AllFreq_neutral$V1, AllFreq_neutral$V2, sep=".")
Loci <- rownames(AllFreq_neutral)
AF2_neutral <- AllFreq[,colnames(AllFreq) %in% Loci]
AllFreq_neutral=AF2_neutral[gsub("\\.", "-", rownames(AF2_neutral)) %in% rownames(Env),]

FilePCANeutral <- paste0(outdir,"/populationStructure/pcaNeutral.rds")
pca_neutral <- rda(AllFreq_neutral[,-1], scale=T)
PCs_neutral <- scores(pca_neutral, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = rownames(AllFreq_neutral), PCs_neutral)

###LOAD THIS FROM FILE
bv <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA_utils/RDAResearchPlan/ordiR2step/best_variables.csv")
best_variables <- bv[,2]
#EnvComp <- Env[best_variables]
EnvComp <- Env[,colnames(Env) %in% best_variables]


print("INTERSECTING COORDINATES AND ENVCOMP AND NEUTRAL PCS")
Variables <- data.frame(Coordinates, EnvComp, PopStruct$PC1, PopStruct$PC2, PopStruct$PC3)

##only components without lat, long
VariablesN <- Variables[,3:ncol(Variables)]

##included
Factors <- colnames(EnvComp)




## FULL MODEL: pRDAfull
formula_string <- paste("AllFreq ~ Lat + Long + PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3 + ", paste(Factors, collapse = " + "))
formula <- as.formula(formula_string)
pRDAfull <- rda(formula, Variables)

RsquareAdj(pRDAfull)$r.squared

#anova(pRDAfull)

##SUBMODEL1: climate model
formula_string <-paste("AllFreq ~ ", paste(Factors, collapse = " + "), "+ Condition(Lat + Long + PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3)")
formula <- as.formula(formula_string)
pRDAclim <- rda(formula, Variables)

RsquareAdj(pRDAclim)$r.squared
#anova(pRDAclim)  

## SUBMODEL 2: population strucutre model
formula_string <- paste("AllFreq ~ PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3 + Condition(Lat + Long +", paste(Factors, collapse = " + "), ")")
formula_struct <- as.formula(formula_string)
pRDAstruct <- rda(formula_struct, Variables)

RsquareAdj(pRDAstruct)$r.squared

 

####################
##Pure geography model
formula_string <- paste("AllFreq ~ Lat + Long + Condition(", paste(Factors, collapse = " + "), "+ PopStruct.PC1 + PopStruct.PC2 + PopStruct.PC3)")
formula_geog <- as.formula(formula_string)
pRDAgeog <- rda(formula_geog, Variables)

#pRDAgeog <- rda(AllFreq ~ Longitude + Latitude + Condition(MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAgeog)$r.squared


##confounded variance
confVar <- RsquareAdj(pRDAfull)$r.squared - RsquareAdj(pRDAclim)$r.squared - RsquareAdj(pRDAstruct)$r.squared - RsquareAdj(pRDAgeog)$r.squared
vp <- varpart(AllFreq, EnvComp, PopStruct[,2:4], Coordinates)
#VennDiagram <- plot(vp)
#ggsave(file = paste0(outdir, "/VariancePartitions",errorval_lower,".png"), VennDiagram, width = 15, height = 6)

png(file = paste0(outdir, "/VarMAFtest_", errorval_lower, ".png"), width = 15, height = 10, units = "in", res = 300)
#plot(vp, Xnames = c("Environment", "PopulationStructure", "Geography"))
plot(vp, bg= c("#01985a", "#009ddb", "#f3c400"), Xnames = c("Environment", "PopulationStructure", "Geography"), digits=3)
dev.off()

result <- rbind(vp$part$fract, vp$part$indfract, vp$part$contr1)
write.csv(result, file=paste0(outdir, "/Var_MAFtest",errorval_lower,".csv"))


###### just tow ork in studio 
library(tidyverse)

# Specify the directory containing the CSV files
#directory <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/RDA"

# List all CSV files in the directory
csv_files <- list.files(path = outdir, full.names = TRUE) # pattern = "\\Var_MAFtest_.*.csv$", 


file_names <- sapply(csv_files, function(x) {
  sub(".*MAF(.*)\\.csv$", "\\1", basename(x))
})

# Read all CSV files into a list of data frames
data_list <- lapply(csv_files, read.csv)
names(data_list) <- file_names

# Create an empty list to hold the R.square values for each entry (i)
r_square_values <- list()

# Iterate over each dataset (j) in the data_list
for (j in seq_along(data_list)) {
  # Extract the R.square column
  r_square_col <- data_list[[j]]$Adj.R.square
  
  # Store the R.square values in the list
  r_square_values[[j]] <- r_square_col
}

#r_square_df <- do.call(cbind, r_square_values)
r_square_df <- data.frame(
  Series = seq_along(r_square_values[[1]]),
  do.call(cbind, r_square_values)
)

colnames(r_square_df) <- c("ID",names(data_list))
rownames(r_square_df) <- data_list[[1]]$X


# Convert the data to long format
r_square_long <- r_square_df %>%
  rownames_to_column(var = "ID") %>%
  pivot_longer(cols = -ID, names_to = "MAF", values_to = "Variance") %>%
  mutate(MAF = as.numeric(MAF)) %>%
  arrange(MAF, ID)

r_square_long <- na.exclude(r_square_long)

plott <- ggplot(r_square_long, aes(x = MAF, y = Variance, color = ID, group = ID)) +
  geom_line() +
  geom_point() +
  labs(title = "Variance Values Across Different MAFs",
       x = "MAF Value",
       y = "Variance Value") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.35)) +  # Use scientific notation if needed
  scale_x_continuous(breaks = pretty(r_square_long$MAF)) +  # Adjust x-axis breaks
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if necessary




#png(file =  "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/RDA/VarPartitions_NonConfounded.png", plott, width = 15, height = 6, units = "in", res = 300)
ggsave(file =  paste0(outdir, "/VarMAFtest_Confounded.png"), plott, width = 15, height = 6)

#r_square_long <- r_square_df %>%
#  as.data.frame() %>%
#  mutate(Index = row_number()) %>%
#  pivot_longer(cols = -Index, names_to = "Dataset", values_to = "R.square")

#r_square_long <- r_square_df %>%
#  pivot_longer(cols = Variance, names_to = "R.square", values_to = "Dataset")
#
## Inspect the long-format data frame
#print(r_square_long)
#
#
#r_square_long <- r_square_long %>%
#  pivot_longer(cols = everything(), names_to = "TimePoint", values_to = "R.square") %>%
#  mutate(Series = factor(rep(1:nrow(r_square_long), each = ncol(r_square_long))))
#
#
## Load necessary library for plotting
#library(ggplot2)
#
## Plot the time series of R.square values
#ggplot(r_square_df, aes(x = Index, y = R.square, color = Dataset)) +
#  geom_line() +
#  labs(title = "Time Series of R.square Values Across Datasets",
#       x = "Index",
#       y = "R.square") +
#  theme_minimal()
####or
## Plot the time series of R.square values
#ggplot(r_square_long, aes(x = Index, y = R.square, color = Dataset)) +
#  geom_line() +
#  labs(title = "Time Series of R.square Values Across Datasets",
#       x = "Index",
#       y = "R.square") +
#  theme_minimal() +
#  scale_y_continuous(labels = scales::scientific) +  # Use scientific notation for small values
#  theme(axis.text.y = element_text(size = 8))  
#