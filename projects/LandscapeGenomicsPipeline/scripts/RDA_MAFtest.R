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
af_file<-args[1]
meta2024<-args[2]
outdir<-args[4]
region<-args[5]
errorperc<- as.numeric(args[6])
dir.create(outdir)

meta <- get_meta_data(meta2024)

af <- read.table(af_file, header=TRUE, na.strings = NaN)
loci <- paste(af$Chr, af$Pos, sep=".")
af_sel <- af %>% dplyr::select(3:ncol(af))
af2 <- as.data.frame(t(af_sel))
colnames(af2) <- loci

AllFreq=af2
rownames(AllFreq) <- gsub("\\.", "-", rownames(AllFreq))
common_elements <- intersect(meta$sampleId, row.names(AllFreq))
# Subset the data frames based on the intersection
matchedMeta  <- meta[meta$sampleId %in% common_elements, ]

errorval_lower <- 1*errorperc/100
errorval_higher <- 1 - errorval_lower

freq_mean <- colMeans(AllFreq)
#AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]
AllFreq <- AllFreq[,-which(freq_mean>=errorval_higher | freq_mean<=errorval_lower)]

## Get Worldclim Data

FileEnv <- paste0(outdir,"/Env.csv")
if (file.exists(FileEnv)) {
  # Load the file
  Env <- readRDS(file=FileEnv)
  message("Env/Biolcim file loaded successfully.")
} else {
  # Generate the file
  Env <- get_worldclim_data(matchedMeta)
  saveRDS(Env, file=FileEnv)
  message("File generated and saved successfully.")
  
}

print("This is Env.")
rownames(Env) <- matchedMeta$sampleId
Env <- na.exclude(Env)
matchedMeta  <- meta[meta$sampleId %in% rownames(Env), ]
Coordinates <- as.data.frame(cbind(matchedMeta$sampleId, as.numeric(matchedMeta$lat), as.numeric(matchedMeta$long)))
colnames(Coordinates) <- c("sampleId", "Latitude", "Longitude")



Env2 <- Env
Env <- scale(Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()


## Recovering scaling coefficients
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')
## Climatic table
Env <- as.data.frame(na.exclude(Env))

new_labels <- function(x) {
  x <- gsub("wc2.1_2.5m_b", "B", x)
  return(x)
}


FileEnvPCA <- paste0(outdir,"/EnvPCA.rds")
if (file.exists(FileEnvPCA)) {
  # Load the file
  Env.PCA <- readRDS(file=FileEnvPCA)
  message("File loaded successfully.")
} else {
  # Generate the file
  Env.PCA <- FactoMineR::PCA(Env)
  saveRDS(Env.PCA, file=FileEnvPCA)
  message("File generated and saved successfully.")
  
}

neutral_SNPS<-args[3]
AllFreq_neutral <- read.table(neutral_SNPS, header = T)
colnames(AllFreq_neutral) <- gsub("\\.", "-", colnames(AllFreq_neutral))

AF2_neutral <- as.data.frame(t(AllFreq_neutral %>% dplyr::select(3:ncol(AllFreq_neutral))))
colnames(AF2_neutral) <- paste(AllFreq_neutral$Chr, AllFreq_neutral$Pos, sep=".")
AllFreq_neutral=AF2_neutral[rownames(AF2_neutral) %in% rownames(Env),]

FilePCANeutral <- paste0(outdir,"/pcaNeutral.rds")

if (file.exists(FilePCANeutral)) {
  # Load the file
  pca_neutral <- readRDS(file=FilePCANeutral)
  message("File loaded successfully.")
} else {
  # Generate the file
  pca_neutral <- rda(AllFreq_neutral[,-1], scale=T) # PCA in vegan uses the rda() call without any predictors
  saveRDS(pca_neutral, file=FilePCANeutral)
  message("File generated and saved successfully.")
  
}


PCs_neutral <- scores(pca_neutral, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = rownames(AllFreq_neutral), PCs_neutral)


####2.Variable selection: forward model building procedure 
Variables <- data.frame(Env)
colnames(Variables) <- new_labels(colnames(Variables))
Parameters <- colnames(Variables)


AllFreq <- AllFreq[rownames(AllFreq) %in% rownames(Env),]

###RDA0 needed??
RDA0 <- rda(AllFreq ~ 1,  Variables) 
formula_string <- paste("AllFreq ~", paste(Parameters, collapse = " + "))
# Convert the string to a formula
formula <- as.formula(formula_string)

# Perform the redundancy analysis (RDA)
RDAfull <- rda(formula, data = Variables)

###PERMUTATION TEST

#####2. FS FUCNTION
#install.packages("packfor", repos = "http://R-Forge.R-project.org")
library(packfor)

global_r2 <- RsquareAdj(RDAfull)$adj.r.squared


RsquareAdj(RDAfull)$adj.r.squared

##determining linear dependencies
##As a rule of the thumb, if √VIF>2 multicollinearity is considered high.


sqrt(vif.cca(RDAfull)) 


#anov2_glob <- anova.cca(SelectedRDA)
#anov2_axis <- anova.cca(SelectedRDA, by="axis")
#anov2_terms <- anova.cca(SelectedRDA, by="terms")

Latitude <- as.numeric(Coordinates$Latitude)
Longitude <- as.numeric(Coordinates$Longitude)
Coords <- cbind(Latitude, Longitude)
PC1 <- PopStruct$PC1
PC2 <- PopStruct$PC2
PC3 <- PopStruct$PC3

##PARTIAL RDA including LAT/LONG and neutral population structure
Vars <- cbind(Variables, Latitude, Longitude, PC1, PC2, PC3 )

##Full model including PCs of neutral genetic strucutre and Coordinate Information 
pRDAfull <- rda(AllFreq ~  . , Vars)
summary(pRDAfull)
RsquareAdj(pRDAfull)$adj.r.squared

#anova(pRDAfull)

##Pure climate model
#pRDAclim <- rda(AllFreq ~  . + Condition(Latitude + Longitude),  Vars[,1:(ncol(Vars)-2)])
pRDAclim <- rda(AllFreq ~  . + Condition(Latitude + Longitude + PC1 + PC2 + PC3), Variables)

summary(pRDAclim)
RsquareAdj(pRDAclim)$adj.r.squared
#anova(pRDAclim)

#############################
##### Pure neutral population structure model  
formula_string <- paste("AllFreq ~ PC1 + PC2 + PC3 + Condition(Longitude + Latitude +", paste(colnames(Variables), collapse = " + "), ")")
# Convert the string to a formula
formula <- as.formula(formula_string)
# Perform the redundancy analysis (RDA)
pRDAstruct <- rda(formula, data = Variables)
#pRDAstruct <- rda(AllFreq ~ PC1 + PC2 + PC3 + Condition(Longitude + Latitude + MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS),  Variables)
RsquareAdj(pRDAstruct)$adj.r.squared
anova(pRDAstruct)


####################
##Pure geography model
formula_string <- paste("AllFreq ~ Longitude + Latitude + Condition(", paste(colnames(Variables), collapse = " + "), " + PC1 + PC2 + PC3)")
# Convert the string to a formula
formula <- as.formula(formula_string)
# Perform the redundancy analysis (RDA)
pRDAgeog <- rda(formula, data = Vars)
#pRDAgeog <- rda(AllFreq ~ Longitude + Latitude + Condition(MAR + EMT + MWMT + CMD + Tave_wt + DD_18 + MAP + Eref + PAS + PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDAgeog)$adj.r.squared


##confounded variance
confVar <- RsquareAdj(pRDAfull)$adj.r.squared - RsquareAdj(pRDAclim)$adj.r.squared - RsquareAdj(pRDAstruct)$adj.r.squared - RsquareAdj(pRDAgeog)$adj.r.squared
vp <- varpart(AllFreq, Variables, PopStruct[,2:4], Coords)
#VennDiagram <- plot(vp)
#ggsave(file = paste0(outdir, "/VariancePartitions",errorval_lower,".png"), VennDiagram, width = 15, height = 6)

png(file = paste0(outdir, "/VarMAFtest_", errorval_lower, ".png"), width = 15, height = 6, units = "in", res = 300)
plot(vp)
dev.off()

result <- rbind(vp$part$fract, vp$part$indfract, vp$part$contr1)
write.csv(result, file=paste0("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/RDA/Var_MAFtest",errorval_lower,".csv"))


###### just tow ork in studio 
library(tidyverse)

# Specify the directory containing the CSV files
directory <- "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/RDA"

# List all CSV files in the directory
csv_files <- list.files(path = directory, pattern = "\\VarMAFtest_.*.csv$", full.names = TRUE)


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

colnames(r_square_df) <- c("Variance",names(data_list))
rownames(r_square_df) <- data_list[[1]]$X


# Convert the data to long format
r_square_long <- r_square_df %>%
  rownames_to_column(var = "ID") %>%
  pivot_longer(cols = -ID, names_to = "MAF", values_to = "Variance") %>%
  mutate(MAF = as.numeric(MAF)) %>%
  arrange(MAF, ID)

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
ggsave(file =  "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData2/results/fullgenome/RDA/VarMAFtest_Confounded.png", plott, width = 15, height = 6)

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