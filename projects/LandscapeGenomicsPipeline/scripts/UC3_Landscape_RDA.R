#install.packages('pegas')
args <- commandArgs(TRUE)

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
library(qvalue)
library(ggpubr)
library(geodata)
library('corrr')
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(gridExtra)



# Functions
get_meta_data <- function(csvfile) {
  meta <- read.csv(csvfile, header = TRUE)
  meta.sub <- meta %>% dplyr::select(sampleId, continent, country, province, lat, long)
  #meta.sub$sampleId <- gsub("\\-", ".", meta.sub$sampleId)
  return(meta.sub)
}

get_worldclim_data <- function(meta.sub) {
  biod <- worldclim_global(var = "bio", 2.5, "data")
  bio <- raster::extract(biod, cbind(meta.sub$long, meta.sub$lat))
  bio.sub <- as.data.frame(bio)
  return(bio.sub)
}

af_file<-args[1]
meta2024<-args[2]
outdir<-args[4]
region<-args[6]

worldclimfolder="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/wc2.5"
envdata=args[5]
dir.create(outdir)

print("LOADING , ALLELE FREQUENCIES AND ANNOTATIONS")
meta <- read.csv(meta2024, header = TRUE)
af <- read.table(af_file, header=TRUE, na.strings = NaN)
#ann="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData/results/annotations.txt"
ann<-args[7]
annotations <- read.table(ann, na.strings = NaN)
annotations <- unique(annotations)
colnames(annotations) <- c("Chr", "Pos", "Gene")
rownames(annotations) <- paste0(annotations$Chr,".", annotations$Pos)
#print(colnames(annotations))
#print(rownames(annotations))
loci <- paste(af$Chr, af$Pos, sep=".")
print(paste("Received", length(loci), "loci.", sep=" "))
##1. Load AF data for the populations of interest 
af2 <- as.data.frame(t(af %>% dplyr::select(3:ncol(af))))
colnames(af2) <- loci


AllFreq=af2
rownames(AllFreq) <- gsub("\\.", "-", rownames(AllFreq))

freq_mean <- colMeans(AllFreq)
AllFreq <- AllFreq[,-which(freq_mean>=0.95 | freq_mean<=0.05)]
print("Working with Allele Frequencies Means >=0.95 or <=0.05 ")
print(nrow(AllFreq))
print(ncol(AllFreq))


print("READING ENVIRONMENTAL DATA")
Env <- read.csv(envdata, header=TRUE)
columns_with_only_na_string <- sapply(Env, function(x) all(x == "na"))
na_string_columns <- names(Env)[columns_with_only_na_string]
Rclean <- Env[, colSums(Env == "na") != nrow(Env)]
Rclean <- as.data.frame(Rclean[2:ncol(Rclean)])
df_numeric <- as.data.frame(lapply(Rclean, as.numeric))

Env2 <- Env

print("READING NULL VALUES DATA")
null_values <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/NUllValuesDetail.csv", header=FALSE)
null_values$V1 <- gsub("--7000--", "..7000..", null_values$V1)



rownames(Env) <- c(Env2$sample)
Env <- as.data.frame(Env)
Rclean <- Env[, colSums(Env == "na") != nrow(Env)]
df_numeric <- as.data.frame(lapply(Rclean, as.numeric))
#Env <- na.exclude(Env)


# Nur Spaltennamen auswählen, die sowohl in null_values$V1 als auch in df_numeric vorhanden sind
valid_columns <- intersect(null_values$V1, colnames(df_numeric))

# Diese validen Spalten aus df_numeric auswählen
new <- df_numeric[, valid_columns]

filtered_null_values <- null_values[null_values$V1 %in% colnames(df_numeric),]
named_list <- as.list(setNames(filtered_null_values$V2, filtered_null_values$V1))

null_values_list <- named_list

df_new <- as.data.frame(lapply(seq_along(new), function(i) {
  replace(new[[i]], new[[i]] == null_values_list[[names(new)[i]]], NA)
}))

colnames(df_new) <- colnames(new)
colnames(df_new) <- gsub("ds.earthserver.xyz..7000..", "", colnames(df_new))
colnames(df_new) <- gsub("global_pesticide_grids_apr", "pest", colnames(df_new))
colnames(df_new) <- gsub("C3S_satellite_soil_moisture", "C3S_sm", colnames(df_new))
colnames(df_new) <- gsub("carbonmonoxide", "", colnames(df_new))
colnames(df_new) <- gsub("carbonmonoxide", "", colnames(df_new))

#na_columns <- names(df_new)[apply(df_new, 2, function(x) all(is.na(x)))]

replace_na_with_median <- function(x) {
  # Compute the median for non-NA values in the row
  row_median <- median(x, na.rm = TRUE)
  # Replace NA values with the median
  x[is.na(x)] <- row_median
  return(x)
}

#df_new_clean <- df_new[, !apply(df_new, 2, function(x) all(is.na(x)))]
rownames(df_new)<- Env2$sample
#core_europe <- df_new[!grepl("^UA|RU|BY|DK", rownames(df_new)), ]


na_count_col <- sapply(df_new, function(x) sum(is.na(x)))
cols_to_keep <- na_count_col <= nrow(df_new)*0.05
remains <- df_new[, cols_to_keep]

print("THESE FACTORS REMAIN")
print(colnames(remains))

# Apply the function to each row
ENV <- as.data.frame(t(apply(remains, 1, replace_na_with_median)))

print("SCALING AND CENTERING")
Env <- scale(ENV, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
scale_env <- attr(Env, 'scaled:scale')
center_env <- attr(Env, 'scaled:center')


print("FIRST PCA")
#colnames(remains) <- paste(substr(colnames(remains), 1, 4), 1:ncol(remains), sep = "_")
pr_eur <- prcomp(Env, scale. = TRUE)
Env.PCA <- princomp(Env)
fviz_pca_var(pr_eur, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

cumvarPlot <- fviz_eig(pr_eur, addlabels = TRUE) + theme_bw()
Env.PCA2 <- pr_eur
## Climatic table
#row.names(Env) <- c(matchedMeta$sampleId)
#Env <- as.data.frame(na.exclude(Env))

new_labels <- function(x) {
  x <- gsub("wc2.1_2.5m_b", "B", x)
  return(x)
}


#Env_PCA <- rda(Env, scale=T)
#screeplot(Env_PCA, type = "barplot", main="PCA Eigenvalues")
#summary_PCA <- summary(Env_PCA)
#cumulative_variance <- summary_PCA$cont$importance[3, ]
#threshold <- 0.90
#num_pcs <- which(cumulative_variance >= threshold)[1]
#FileEnvPCA <- paste0(outdir,"/EnvPCA.rds")
#if (file.exists(FileEnvPCA)) {
#  # Load the file
#  Env.PCA <- readRDS(file=FileEnvPCA)
#  message("File loaded successfully.")
#} else {
#  # Generate the file
#  Env.PCA <- princomp(Env)
#  saveRDS(Env.PCA, file=FileEnvPCA)
#  message("File generated and saved successfully.")
#}

#Env.PCA <- princomp(Env)
cumvarPlot <- fviz_eig(Env.PCA, addlabels = TRUE) + theme_bw()
Env.PCA2 <- Env.PCA
#rownames(Env.PCA2$loadings) <- new_labels(rownames(Env.PCA2$loadings))
RDAPlot <- fviz_pca_var(Env.PCA2, col.var = "black") + theme_bw() 
ggsave(paste0(outdir,"/CumVar.png"), plot = cumvarPlot, width = 8, height = 6, dpi = 300)
ggsave(paste0(outdir,"/EnvRDA.png"), plot = RDAPlot, width = 8, height = 6, dpi = 300)


###In PCA, cos2 refers to the squared cosine of the variables on the principal components.
#It represents the quality of representation of each variable on the principal components.
#Higher values of cos2 indicate that the variable is well-represented by the principal component.


# Create an empty list to store the plots
plot_list <- list()

# Loop through axes 1 to 4 and store each plot
for (i in 1:4) {
  plot_list[[i]] <- fviz_cos2(Env.PCA2, choice = "var", axes = i) +
    scale_x_discrete() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)  # Rotate x-axis labels
    )
}

# Display the plots in a grid (2 rows and 2 columns)
#grid.arrange(grobs = plot_list, nrow = 2, ncol = 2)

# Save the grid to a PNG file
ggsave(paste0(outdir,"/VariablePartitions.png"), width = 8, height = 6, dpi = 300)
# Use grid.arrange inside grid.draw to save it
grid.arrange(grobs = plot_list, nrow = 2, ncol = 2)
# Turn off the device to finalize the file
dev.off()

x <- fviz_pca_var(Env.PCA2, col.var = "cos2", gradient.cols = c("black", "orange", "green"), repel = TRUE) + theme_bw()
# Save the plot
ggsave(paste0(outdir, "/EnvRDA_Colored.png"), plot = x, width = 8, height = 6, dpi = 300)


summary_pca <- summary(pr_eur)
write.table(capture.output(summary_pca), file = paste0(outdir,"/SummaryPCA0.txt"), sep = "\t", quote = FALSE)

explained_variance <- summary_pca$sdev^2 / sum(summary_pca$sdev^2)
cumulative_variance <- cumsum(explained_variance)
threshold <- 0.80  
num_pcs <- which(cumulative_variance >= threshold)[1]


#EnvComp <- as.data.frame(pr_eur$scores[,1:num_pcs])
EnvComp <- as.data.frame(Env.PCA$scores[,1:num_pcs])
#Coordinates <- metadata[,1:3]
#Coordinates <- Coordinates[Coordinates$sample %in% Env2$sample,]
Coordinates <- as.data.frame(cbind(meta$sampleId, meta$lat, meta$long))
colnames(Coordinates) <- c("sample", "lat", "long")
#print(Coordinates$sample)
#print(Env2$sample)
Coordinates <- Coordinates[Coordinates$sample  %in% rownames(EnvComp) ,]
#print(Coordinates)

print("ANALYSING ENVIRONMENTAL DATA COMPOSITION")
print(nrow(EnvComp))
print(nrow(Coordinates))
EnvComp <- EnvComp[rownames(EnvComp) %in% Coordinates$sample, ]

#print(setdiff(Coordinates$sample, rownames(EnvComp))) 
#print(setdiff(rownames(EnvComp), Coordinates$sample)) 

Env.PCA <- pr_eur
cos2_values <- get_pca_var(Env.PCA)$cos2
best_vars <- apply(cos2_values, 2, function(x) rownames(cos2_values)[which.max(x)])
topvars <- as.character(best_vars[1:num_pcs])

###ADD PHENOTYPIC DATA?

print("PREPARING RDA NEUTRAL")
###Population Structure: PCA on intergenic SNPs
neutral_SNPS<-args[3]
#neutral_SNPS="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/FullData/results/fullgenome2/Neutral.final.af"

AllFreq_neutral <- read.table(neutral_SNPS, header = T)
AF2_neutral <- as.data.frame(t(AllFreq_neutral %>% dplyr::select(3:ncol(AllFreq_neutral))))
colnames(AF2_neutral) <- paste(AllFreq_neutral$Chr, AllFreq_neutral$Pos, sep=".")
#AF3_neutral <- AF2_neutral[3:nrow(AF2_neutral),]
#af10k <- as.data.frame(af3[,1:10000])

rownames(Env) <- c(Env2$sample)

#AllFreq=af10k
AllFreq_neutral=AF2_neutral[gsub("\\.", "-", rownames(AF2_neutral)) %in% rownames(Env),]

print("PERFORMING RDA NEUTRAL")

## Running a PCA on neutral genetic markers
#pca_neutral <- rda(AllFreq_neutral[,-1], scale=T) # PCA in vegan uses the rda() call without any predictors
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

png(filename = paste0(outdir,"/RDA_PopStructure.png"), width = 800, height = 600)
screeplot(pca_neutral, type="barplot")
#plot(pca_neutral)
dev.off()

PCs_neutral <- scores(pca_neutral, choices=c(1:3), display="sites", scaling=0)

PopStruct <- data.frame(Population = rownames(AllFreq_neutral), PCs_neutral)
#colnames(PopStruct) <- c("Population", "PC1", "PC2", "PC3")#
###These PCs need to be added to "Variables" in order to get outlier loci later


print("INTERSECTING COORDINATES AND ENVCOMP")
###No variable selection needed?! becasue environemtnal PCAs
### Variable Selection
Variables <- data.frame(Coordinates, EnvComp)
###reduce to AF rows
Variables <- Variables[Variables$sample %in% rownames(AllFreq),]
#admin <- ne_countries(scale = "medium", returnclass = "sf")

#Variables$Latitude <- as.numeric(Variables$Latitude)
#Variables$Longitude <- as.numeric(Variables$Longitude)
VariablesN <- Variables[,4:ncol(Variables)]

##included
column_names <- colnames(VariablesN)
print(colnames(VariablesN))
# Filter column names to include only those starting with "Comp."
comp_columns <- grep("^Comp\\.", column_names, value = TRUE)
# Create the formula string dynamically
formula_string <- paste("AllFreq ~ lat + long +", paste(comp_columns, collapse = " + "))
print(formula_string)
# Convert the string to a formula
formula <- as.formula(formula_string)
AllFreq <- AllFreq[rownames(AllFreq) %in% Variables$sample,]
print(nrow(AllFreq))


print("PERFORMING RDA FULL")
#Variance Partitioning 
#pRDAfull <- rda(AllFreq ~ Variables$Latitude + Variables$Longitude + Variables$Comp.1 + Variables$Comp.2 + Variables$Comp.3 + Variables$Comp.4,  Variables)
#pRDAfull <- rda(AllFreq ~ Latitude + Longitude + Comp.1 + Comp.2 + Comp.3,  Variables)
# Run the rda function with the dynamically created formula
FilepRDAfull <- paste0(outdir,"/pRDAfull.rds")
if (file.exists(FilepRDAfull)) {
  # Load the file
  pRDAfull <- readRDS(file=FilepRDAfull)
  message("File 'pRDAfull.rds' loaded successfully.")
} else {
  # Generate the file
  pRDAfull <- rda(formula, Variables)
  saveRDS(pRDAfull, file=FilepRDAfull)
  message("File generated and saved successfully.")
  
}

#pRDAfull <- rda(formula, Variables)
RsquareAdj(pRDAfull)
#anova(pRDAfull)
#save(pRDAfull, file=paste0(outdir,"/pRDAfull.rds"))

png(filename = paste0(outdir,"/pRDA_full.png"), width = 800, height = 600)
plot(pRDAfull)
dev.off()

print("PERFORMIGN RDA: CLIMATE MODEL")
##Pure climate model

#pRDAclim <- rda(AllFreq ~ Comp.1 + Comp.2 + Comp.3 + Comp.4 + Condition(Latitude + Longitude),  Variables)
formula_string <- paste("AllFreq ~ ", paste(comp_columns, collapse = " + "))
formula_string <-paste("AllFreq ~ ", paste(comp_columns, collapse = " + "), "+ Condition(lat + long)")
#pRDAclim <- rda(AllFreq ~ Comp.1 + Comp.2 + Comp.3 + Comp.4, VariablesN)
formula <- as.formula(formula_string)

FilepRDAclim <- paste0(outdir,"/pRDAclim.rds")
if (file.exists(FilepRDAclim)) {
  # Load the file
  pRDAclim <- readRDS(file=FilepRDAclim)
  message("File loaded successfully.")
} else {
  # Generate the file
  pRDAclim <- rda(formula, Variables)
  saveRDS(pRDAclim, file=FilepRDAclim)
  message("File generated and saved successfully.")
}

#pRDAclim <- rda(formula, Variables)
RsquareAdj(pRDAclim)
#anova(pRDAclim)

png(filename = paste0(outdir,"/pRDA_climate.png"), width = 800, height = 600)
plot(pRDAclim)
dev.off()


png(filename = paste0(outdir,"/CorrPlot.png"), width = 800, height = 600)
corrplot(cor(VariablesN), type="upper")
dev.off()

#Genotype-Environment Associations: identifying loci under selection with Neutral Structure as condition
## NOT DONE YET
#formula_string <- paste("AllFreq ~ ", paste(comp_columns, collapse = " + "), "Condition")

PC_names <- colnames(PCs_neutral)
formula_string <-paste("AllFreq ~ ", paste(comp_columns, collapse = " + "), "+ Condition(", paste(PC_names, collapse = " + "), ")")
formula_env <- as.formula(formula_string)
Var2 <- cbind(Variables, PCs_neutral)
print("Making RDAEnv")

FileRDA_ENV <- paste0(outdir,"/RDA_env.rds")
if (file.exists(FileRDA_ENV)) {
  # Load the file
  RDA_env <- readRDS(file=FileRDA_ENV)
  message("File loaded successfully.")
} else {
  # Generate the file
  RDA_env <- rda(formula_env, Var2)
  saveRDS(RDA_env, FileRDA_ENV)
  message("File generated and saved successfully.")
}

#RDA_env <- rda(formula_env, Var2)
#placehodlderRDA RDA_env <- rda(AllFreq ~ Comp.1 + Comp.2 + Comp.3 + Comp.4, VariablesN)
#screeplot(RDA_env, main="Eigenvalues of constrained axes")

source("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/src/rdadapt.R")

###conduct rdadapt_env
rdadapt_env<-rdadapt(RDA_env, 2)

## P-values threshold after Bonferroni correction
thres_env <- 0.1/length(rdadapt_env$p.values)

png(filename = paste0(outdir,"/Hist_rdadapt_pvalues.png"), width = 800, height = 600)
hist(rdadapt_env$p.values)
dev.off()

## Identifying the loci that are below the p-value threshold
outliers <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], p.value = rdadapt_env$p.values[which(rdadapt_env$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env$p.values<thres_env)], split = "\\."), function(x) x[1])))

print("THESE ARE THE OUTLIERS:")
print(nrow(outliers))
#print(outliers)


## Top hit outlier per contig
outliers <- outliers[order(outliers$contig, outliers$p.value),]
## List of outlier names, one outlier per contig/chromosome
##outliers_rdadapt_env <- as.character(outliers$Loci[!duplicated(outliers$contig)])

outliers_rdadapt_env <- as.character(outliers$Loci)

#top_outliers_per_chromosome <- outliers %>%
#  group_by(contig) %>%
#  slice_head(n = 10) %>%
#  ungroup()
top_1perc_outliers_per_chromosome <- outliers %>%
  group_by(contig) %>%
  filter(rank(p.value) <= round(n() * 0.1)) %>% # Selecting top 10% outliers per chromosome
  ungroup()


#top_chromosome_outliers <- as.character(top_outliers$Loci)
top_1p_chromosome_outliers <- as.character(top_1perc_outliers_per_chromosome$Loci)
print("Top 1% outleirs:")
print(nrow(top_1perc_outliers_per_chromosome))

## Formatting table for ggplot
locus_scores <- scores(RDA_env, choices=c(1:2), display="species", scaling="none") # vegan references "species", here these are the loci
TAB_loci <- data.frame(names = rownames(locus_scores), locus_scores)
TAB_loci$chromosome <- sub("\\..*", "", TAB_loci$names)
TAB_loci$type <- "Neutral"
TAB_loci$type[TAB_loci$names%in%outliers$Loci] <- "All outliers"
#TAB_loci$type[TAB_loci$names%in%outliers_rdadapt_env] <- "Top outliers"
##TAB_loci$type[TAB_loci$names%in%top_chromosome_outliers] <- "Top outliers"
TAB_loci$type[TAB_loci$names%in%top_1p_chromosome_outliers] <- "Top outliers"
TAB_loci$type <- ifelse(TAB_loci$type == "Top outliers", 
                        paste0("Chr ", TAB_loci$chromosome), 
                        TAB_loci$type)
#TAB_loci$type <- factor(TAB_loci$type, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_loci$type <- factor(TAB_loci$type, levels = c(unique(TAB_loci$type)))
#TAB_loci["Gene"] <- annotations$Gene
TAB_loci <- TAB_loci[order(TAB_loci$type),]
TAB_var <- as.data.frame(scores(RDA_env, choices=c(1,2), display="bp")) # pull the biplot scores

#print(TAB_loci)
TAB_loci$Gene <- annotations$Gene[rownames(annotations) %in% rownames(TAB_loci)]

## Biplot of RDA loci and variables scores
PLOT1 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*20, y=RDA2*20, colour = type), size = 1.4) +
  geom_text_repel(data = TAB_loci %>% filter(names %in% top_1p_chromosome_outliers), aes(x=RDA1*20, y=RDA2*20, label=Gene), size=2.5, vjust=1.5, hjust=0.5, family="Times") + # Label only top outliers
  ##scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  scale_color_manual(values=c( "gray90", "#F9A242FF", "#44ad61","#7761d0","#4ab09c","#967dca","#588dcc")) + 
  #scale_color_manual(values=c( "gray90", "#F9A242FF", "#44ad61")) +
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab("RDA 1") + ylab("RDA 2") +
  facet_wrap(~"RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") 
  #theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

ggsave(paste0(outdir,"/RDA_comp.png"), plot = PLOT1, width = 8, height = 6, dpi = 300)

## Manhattan plot
Outliers <- rep("Neutral", length(colnames(AllFreq)))
Outliers[colnames(AllFreq)%in%outliers$Loci] <- "All outliers"
Outliers[colnames(AllFreq)%in%outliers_rdadapt_env] <- "Top outliers"
Outliers <- factor(Outliers, levels = c("Neutral", "All outliers", "Top outliers"))
TAB_manhatan <- data.frame(pos = 1:length(colnames(AllFreq)), 
                           pvalues = rdadapt_env$p.values, 
                           Outliers = Outliers,
                           loci = colnames(AllFreq)) %>% separate(loci, into = c('Chr', 'Pos'), sep = '\\.')
#TAB_manhatan <- TAB_manhatan %>% separate(loci, into = c('Chr', 'Pos'), sep = '\\.')
TAB_manhatan <- TAB_manhatan[order(TAB_manhatan$Outliers),]
PLOT1M <- ggplot(data = TAB_manhatan) +
  geom_point(aes(x=pos, y=-log10(pvalues), col = Outliers), size=1.4) +
  #geom_text_repel(data = TAB_manhatan %>% filter(names %in% top_1p_chromosome_outliers), aes(x=RDA1*20, y=RDA2*20, label=Gene), size=2.5, vjust=1.5, hjust=0.5, family="Times") + # Label only top outliers
  facet_grid(. ~ Chr, scales = "free_x", space = "free") + 
  scale_color_manual(values = c("gray90", "#F9A242FF", "#6B4596FF")) +
  xlab("Loci") + ylab("-log10(p.values)") +
  geom_hline(yintercept=-log10(thres_env), linetype="dashed", color = gray(.80), size=0.6) +
  #facet_wrap(~"Manhattan plot", nrow = 3) +
  ggtitle("Manhattan plot") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") 
  #theme(legend.position="right", legend.background = element_blank(), panel.grid = element_blank(), legend.box.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11))

ggsave(paste0(outdir,"/RDA_comp_Manhattan_split.png"), plot = PLOT1M, width = 8, height = 6, dpi = 300)


#compare outliers

#not account for pop structure
## Running a simple RDA model

formula_string <- paste("AllFreq ~ ", paste(comp_columns, collapse = " + "))
formula <- as.formula(formula_string)
RDA_env_unconstrained <- rda(formula, Variables)
#RDA_env_unconstrained <- rda(AllFreq ~ Comp.1 + Comp.2 + Comp.3 ,  Variables)

rdadapt_env_unconstrained <- rdadapt(RDA_env_unconstrained, 2)
thres_env <- 0.1/length(rdadapt_env_unconstrained$p.values)

## Identifying the outliers for the simple RDA
outliers_unconstrained <- data.frame(Loci = colnames(AllFreq)[which(rdadapt_env_unconstrained$p.values<thres_env)], p.value = rdadapt_env_unconstrained$p.values[which(rdadapt_env_unconstrained$p.values<thres_env)], contig = unlist(lapply(strsplit(colnames(AllFreq)[which(rdadapt_env_unconstrained$p.values<thres_env)], split = "_"), function(x) x[1])))
outliers_unconstrained <- outliers_unconstrained[order(outliers_unconstrained$contig, outliers_unconstrained$p.value),]
outliers_rdadapt_env_unconstrained <- as.character(outliers_unconstrained$Loci[!duplicated(outliers_unconstrained$contig)])

## For all the outliers
list_outliers_RDA_all <- list(RDA_constrained = as.character(outliers$Loci), RDA_unconstrained = as.character(outliers_unconstrained$Loci))
ggVennDiagram(list_outliers_RDA_all, category.names = c("partial RDA", "simple RDA"), lty="solid", size=0.2) + 
  scale_fill_gradient2(low = "white", high = 'gray40') + scale_color_manual(values = c("grey", "grey", "grey", "grey")) + guides(fill = "none") + theme(text = element_text(size=16, family = "Times"))


##
##
##
list_outliers_RDA_top <- list(RDA_constrained = outliers_rdadapt_env, RDA_unconstrained = outliers_rdadapt_env_unconstrained)
common_outliers_RDA_top <- Reduce(intersect, list_outliers_RDA_top)

##adaptively enriched genetic space
formula_string <- paste("AllFreq[,common_outliers_RDA_top] ~ ", paste(comp_columns, collapse = " + "))
formula <- as.formula(formula_string)
RDA_outliers <- rda(formula, Variables)
#RDA_outliers <- rda(AllFreq[,common_outliers_RDA_top] ~ Comp.1 + Comp.2 + Comp.3 ,  Variables)

## RDA biplot ###NOT WORKING YET
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
ex1<- as.character(round(RDA_outliers$CCA$eig[1],1))
ex2<- as.character(round(RDA_outliers$CCA$eig[2],1))
PLOT2 <- ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#EB8055FF", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", size=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "Times") +
  xlab(paste0("RDA 1 (", ex1,")%")) + ylab(paste0("RDA 2 (", ex2,")%")) +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "Times") 
  #theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))

ggsave(paste0(outdir,"/RDA_outliers.png"), plot = PLOT2, width = 8, height = 6, dpi = 300)

## Loading the climatic rasters




##### ENVIRONMETNAL STACK
#
#
#
#ras <- stack(list.files(worldclimfolder, pattern = ".tif$", full.names = T))
##names(ras) <- unlist(strsplit(unlist(lapply(strsplit(list.files("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/wc2.1/", pattern = ".tif$", full.names = T), split = "/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/wc2.1/wc2.1_10m_"), function(x) x[2])), split = ".tif"))
##ras_6190 <- ras[[grep("6190", names(ras))]]
###
##bio1 <- ras[[1]]
##bio4 <- ras[[14]]
##bio7 <- ras[[17]]
#
#selected_layers_list <- list()
#
## Loop through each value in best_var_value
#for (value in bioclims) {
#  # Extract layers with names matching the current value
#  selected_layers_list[[value]] <- ras[[value]]
#}
#
#Env_Stack_wc <- stack(selected_layers_list)
#
#
#remove.NAs.stack<-function(rast.stack){
#  nom<-names(rast.stack)
#  test1<-calc(rast.stack, fun=sum)
#  test1[!is.na(test1)]<-1
#  test2<-rast.stack*test1
#  test2<-stack(test2)
#  names(test2)<-nom
#  return(test2)
#}
#
##bio1 <- remove.NAs.stack(bio1)
##desired_variable <- "Near.Surface.Air.Temperature"  # Update this with your desired variable name
##extracted_data <- nc_file[[desired_variable]]
#
##nc <- raster(nc_file)
##plot(nc, main = "Near Surface Air Temperature", col = rev(terrain.colors(100)))
##Env_Stack <- stack(bio1, bio12, bio4, bio7)
#Env_Stack_polished <- remove.NAs.stack(Env_Stack_wc)
#europe_extent <- extent(-25, 45, 34, 72)
##north_america_extent <- extent(-172, -10, 7, 83)
#
## Crop the raster to the extent of Europe
#europe_data <- crop(Env_Stack_polished, europe_extent)
##NA_data <- crop(nc, north_america_extent)
##plot(europe_data, main = "Envs", col = rev(terrain.colors(100)))
##plot(NA_data, main = "Near Surface Air Temperature in North America", col = rev(terrain.colors(100)))
#
#world <- ne_countries(scale = "medium", returnclass = "sf")
#europe <- subset(world, continent == region)
##WE <- subset(world, subregion == "Western Europe")
##Austria <- subset(world, name_en == "Austria")
#
#Coordinates <- cbind(Variables$Latitude, Variables$Longitude)
##Env <- cbind(meta2$NSAT, meta2$BIO1, meta2$nFlies)
#
###new
#selected_columns_list <- list()
#
## Loop through each value in best_var_value
#for (value in bioclims) {
#  # Extract columns with names matching the current value
#  selected_columns_list[[value]] <- Env2[, value, drop = FALSE]
#}
#
#Env_Stack_Env <- do.call(cbind, selected_columns_list)
#
#
##Env_Stack_Env <- cbind(Env$wc2.1_2.5m_bio_1, Env$wc2.1_2.5m_bio_4, Env$wc2.1_2.5m_bio_7, Env$wc2.1_2.5m_bio_12)
##colnames(Env_Stack_Env) <-  names(Env_Stack_polished)
#
#Env <- scale(Env_Stack_Env, center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()
#scale_env <- attr(Env, 'scaled:scale')
#center_env <- attr(Env, 'scaled:center')
#
#row.names(RDA_outliers$CCA$biplot) <- names(Env_Stack_polished)
#
#admin <- ne_countries(scale = "medium", returnclass = "sf")
#range<-europe
####adaptive index- WORKS IF ENV NAMES ARE ENV NAMES from RDA and scale names adn so on match
#source("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/src/adaptive_index.R")
#res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = Env_Stack_polished, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
##res_RDA_proj_WE <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = Env_Stack_polished, range = WE, method = "loadings", scale_env = scale_env, center_env = center_env)
##res_RDA_proj_AUT <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = Env_Stack_polished, range = Austria, method = "loadings", scale_env = scale_env, center_env = center_env)
#
####excursion
#
#
#
#
#
####BEFORE
##
##range_polygons <- Polygons(list(range_polygon), "range")
##range_spatial_polygons <- SpatialPolygons(list(range_polygons))
##range <- SpatialPolygonsDataFrame(range_spatial_polygons, data.frame(id = 1), match.ID = FALSE)
##crs(range) <- CRS("+proj=longlat +datum=WGS84")
##writeOGR(range, dsn = "./Data/pinucon", layer = "pinucont", driver = "ESRI Shapefile")
##res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 1, env_pres = NSAT, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
###### BEFORE END
#
#
#
### Vectorization of the climatic rasters for ggplot
#RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
##RDA_proj <- list(res_RDA_proj_WE$RDA1, res_RDA_proj_WE$RDA2)
##RDA_proj <- list(res_RDA_proj_AUT$RDA1, res_RDA_proj_AUT$RDA2)
#RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
#
#for(i in 1:length(RDA_proj)){
#  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
#}
#
#
#europe <- admin[admin$continent == "Europe", ]
#europe_simplified <- rmapshaper::ms_simplify(europe, keep = 0.01)
#europe_simplified <- europe_simplified[, c("name", "geometry")]
#
####
### Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
#TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
#colnames(TAB_RDA)[3] <- "value"
#TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
#AD_INDEX <- ggplot(data = TAB_RDA) + 
#  #geom_sf(data = admin, fill=gray(.9), size=0) +
#  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) +
#  coord_cartesian(xlim = c(-10, 40), ylim = c(35,60)) +
#  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
#  #geom_sf(data = europe_simplified, fill=NA, size=0.1) +
#  #coord_sf(xlim = c(-148, -98), ylim = c(35, 64), expand = FALSE) +
#  xlab("Longitude") + ylab("Latitude") +
#  guides(fill=guide_legend(title="Adaptive index")) +
#  facet_grid(~ variable) +
#  theme_bw(base_size = 11, base_family = "Times") +
#  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
#
#ggsave(paste0(outdir,"/AD_Index2.png"), plot = AD_INDEX, width = 12, height = 6, dpi = 300)
#
#### Predicting future genomic offset 
###we need projected data for that 
#
##source("/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/RDA-landscape-genomics/src/genomic_offset.R")
###DATA NEEDED
##res_RDA_proj2080 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_6190, env_fut = ras_2080, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
##res_RDA_proj2050 <- genomic_offset(RDA_outliers, K = 2, env_pres = ras_6190, env_fut = ras_2050, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)
##
#### Table global genetic offset predicted for 2050 and 2080
##RDA_proj_offset <- data.frame(rbind(rasterToPoints(res_RDA_proj2050$Proj_offset_global), rasterToPoints(res_RDA_proj2080$Proj_offset_global)), Date = c(rep("2050", nrow(rasterToPoints(res_RDA_proj2050$Proj_offset_global))), rep("2080", nrow(rasterToPoints(res_RDA_proj2080$Proj_offset_global)))))
##
#### Projecting genomic offset on a map
##colors <- c(colorRampPalette(brewer.pal(11, "Spectral")[6:5])(2), colorRampPalette(brewer.pal(11, "Spectral")[4:3])(2), colorRampPalette(brewer.pal(11, "Spectral")[2:1])(3))
##ggplot(data = RDA_proj_offset) + 
##  geom_sf(data = admin, fill=gray(.9), size=0) +
##  geom_raster(aes(x = x, y = y, fill = cut(Global_offset, breaks=seq(1, 8, by = 1), include.lowest = T)), alpha = 1) + 
##  scale_fill_manual(values = colors, labels = c("1-2","2-3","3-4","4-5","5-6","6-7","7-8"), guide = guide_legend(title="Genomic offset", title.position = "top", title.hjust = 0.5, ncol = 1, label.position="right"), na.translate = F) +
##  geom_sf(data = admin, fill=NA, size=0.1) +
##  coord_sf(xlim = c(-148, -98), ylim = c(35, 64), expand = FALSE) +
##  xlab("Longitude") + ylab("Latitude") +
##  facet_grid(~ Date) +
##  theme_bw(base_size = 11, base_family = "Times") +
##  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
##
#
#
#
####plot PCAs
## Function to read frequency data
#read_frequency_data <- function(filepath) {
#  ## read the allele frequency dataset
#  DATA <- read.table(filepath, header = TRUE, comment.char = "")
#  return(DATA)
#}
#
## Function to process inversion data
#process_inversion_data <- function(DATA, meta.sub, Chr, Start, End) {
#  ## transpose the allele frequeny matrix and only retain SNPs, within the inverted genomic region.
#  DATA.inv <- as.data.frame(t(DATA %>%
#                                filter(DATA$Chr == DATA$Chr & DATA$Pos > Start & DATA$Pos < End) %>%
#                                dplyr::select(3:ncol(DATA))))
#  print(rownames(DATA.inv))
#  ## make new column from rownames (sampleIds)
#  DATA.inv$sampleId <- rownames(DATA.inv)
#  return(DATA.inv)
#}
#
#INV.full <- c("IN2Lt", "IN3RP")
#Chr.full <- c("2L", "3R")
#Start.full <- c(2225744,16432209)
#End.full <- c(13154180,24744010)
#
#
## Function to process non-inversion data
#process_non_data <- function(DATA, meta.sub, Chr.full, Start.full, End.full) {
#  DATA.non.1 <- DATA %>% filter(!((DATA$Chr == Chr.full[1] & DATA$Pos > Start.full[1] & DATA$Pos < End.full[1])))
#  DATA.non <- as.data.frame(t(DATA.non.1 %>% filter(!(DATA.non.1$Chr == Chr.full[2] & DATA.non.1$Pos > Start.full[2] & DATA.non.1$Pos < End.full[2])) %>%
#                                dplyr::select(3:ncol(DATA.non.1))))
#  DATA.non$sampleId <- rownames(DATA.non)
#  return(DATA.non)
#}
#
#
#perform_pca <- function(DATA) {
#  PCA.result <- PCA(DATA, graph = FALSE)
#  return(PCA.result)
#}
#
#map_levels_to_numbers <- function(vec, start, end) {
#unique_levels <- unique(vec)
#level_mapping <- setNames(rep(start:end, length.out = length(unique_levels)), unique_levels)
#num_vec <- level_mapping[vec]
#return(num_vec)
#}
#
## Function to create PCA plot
#create_pca_plot <- function(PCA.result, DATA, region, INV, inside = TRUE) {
#  plot_type <- ifelse(inside, "inside", "outside")
#  if (region == "Europe") {
#    COLOR <- DATA$country
#  } else {
#    COLOR <- DATA$province
#  }
#  if (region == "NorthAmerica") {
#    REG <- "NorthAm."
#  } else {
#    REG <- region
#  }
#  COLOR2 <- map_levels_to_numbers(COLOR, 1, 11)
#  SHAPE <- map_levels_to_numbers(COLOR, 15, 20)
#  PLOT <- fviz_pca_ind(PCA.result,
#                       col.ind = COLOR,
#                       pointsize = 3,
#                       alpha = 0.7,
#                       mean.point = FALSE,
#                       label = COLOR,
#                       repel = TRUE
#  ) +
#    theme_bw() +
#    ggtitle(paste0("PCA - ", REG, " ", plot_type, " ", INV)) +
#    labs(color = "Count(r)y") +
#    labs(shape = "Count(r)y") +
#    guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
#    guides(shape = guide_legend(nrow = 3, byrow = TRUE)) +
#    scale_shape_manual(values = SHAPE) +
#    scale_color_manual(values = COLOR2)
#  return(PLOT)
#}
#
#####TESTING INVERSIONS####
###################
#
#
##FULL <- af
##rownames(FULL) <- paste0(FULL$Chr,".", FULL$Pos)
##random_SNPs <- sample(rownames(FULL), 2700)
##FULL <- as.data.frame(t(FULL[,3:ncol(FULL)]))
### Subset the data frame using the selected column indices
##FULL_subset <- FULL[,random_SNPs]
##FULL.NULL <- PCA(FULL_subset)
##FULL.NULL.ALL <- PCA(FULL)
##FNULL <- create_pca_plot(FULL.NULL, matchedMeta, region, " Random 2.7k", inside = TRUE)
##FULL.PLOT <- create_pca_plot(FULL.NULL.ALL, matchedMeta, region, "All SNPs", inside = TRUE)
##ggsave(file = paste0(outdir, "/Random2kSNPS_PCA.png"), FNULL, width = 15, height = 6)
##ggsave(file = paste0(outdir, "/FullModel_PCA.png"), FULL.PLOT, width = 15, height = 6)
##
##
##af_file_neutral=args[3]
##NEUTRAL <- read_frequency_data(af_file_neutral)
##rownames(NEUTRAL) <- paste0(NEUTRAL$Chr,".", NEUTRAL$Pos)
##
##NEUTRAL.N <- as.data.frame(t(NEUTRAL[,3:ncol(NEUTRAL)]))
##NEUTRAL.NULL <-PCA(NEUTRAL.N)
##NNN <- create_pca_plot(NEUTRAL.NULL, matchedMeta, region, "none", inside = TRUE)
##ggsave(file = paste0(outdir, "/NeutralStruct_NULL.png"), NNN, width = 15, height = 6)
##
##NEUTRAL.inv.IN2Lt <- process_inversion_data(NEUTRAL, matchedMeta, Chr.full[1], Start.full[1], End.full[1])
##PCA.NEUTRAL.inv.IN2Lt <- perform_pca(NEUTRAL.inv.IN2Lt[, 1:(ncol(NEUTRAL.inv.IN2Lt) - 3)])
##PLOT.inv.IN2Lt <- create_pca_plot(PCA.NEUTRAL.inv.IN2Lt, matchedMeta, region, "IN2Lt", inside = TRUE)
##ggsave(file = paste0(outdir,"/NeutralStruct_Inv2Lt.png"), PLOT.inv.IN2Lt, width = 15, height = 6)
##
##NEUTRAL.inv.IN3RP <- process_inversion_data(NEUTRAL, matchedMeta, Chr.full[2], Start.full[2], End.full[2])
##PCA.inv.IN3RP <- perform_pca(NEUTRAL.inv.IN3RP[, 1:(ncol(NEUTRAL.inv.IN3RP) - 3)])
###save_pca_results(PCA.inv.IN3RP, DATA.inv.IN3RP, region, "IN3RP", inside = TRUE)
##PLOT.inv.IN3RP <- create_pca_plot(PCA.inv.IN3RP, matchedMeta, region, "IN3RP", inside = TRUE)
##ggsave(file = paste0(outdir,"/NeutralStruct_Inv3RP.png"), PLOT.inv.IN3RP, width = 15, height = 6)
##
##NEUTRAL.non <- process_non_data(NEUTRAL, matchedMeta, Chr.full, Start.full, End.full)
##PCA.non <- perform_pca(NEUTRAL.non[, 1:(ncol(NEUTRAL.non) - 3)])
##PLOT.non <- create_pca_plot(PCA.non, matchedMeta, region, "NonInversions", inside = FALSE)
##ggsave(file = paste0(outdir,"/NeutralStruct_nonInversion.png"), PLOT.non, width = 15, height = 6)
##
##for (i in seq(1, ncol(FULL.NULL$var$contrib))) {
##  # Create data frame for the current dimension
##  Dim1 <- as.data.frame(FULL.NULL$var$contrib[, i])
##  # Add Locus as a column
##  Dim1$Locus <- rownames(Dim1)
##  # Separate Locus into Chr and Pos
##  DD <- Dim1 %>% separate(Locus, into = c('Chr', 'Pos'), sep = '\\.')
##  # Rename columns appropriately
##  colnames(DD)[1] <- "P"
##  # Create plot
##  PLOT <- ggplot(DD, aes(x = Pos, y = P)) +
##    geom_point(col = rgb(0, 0, 0, 0.4), pch = 16) +
##    facet_grid(. ~ Chr, scales = "free_x", space = "free") +
##    xlab("Genomic Position [Mbp]") +
##    ylab(paste0("Contribution Dim", i)) +
##    theme_bw()
##  # Save the plot
##  ggsave(file = paste0(outdir, "/PopStructPCA_contribution_Dim", i, ".png"), PLOT, width = 15, height = 6)
##}
#
#
########## outlier check for both RDA 1 adn 2 but focus in 2
#
#lociscores <- as.data.frame(scores(RDA_outliers, display="species"))
#lociscores['loci'] <- rownames(lociscores)
#mean_score_RDA1 <- mean(lociscores$RDA1)
#sd_score_RDA1 <- sd(lociscores$RDA1)
#threshold_RDA1 <- mean_score + 3 * sd_score
#threshold_RDA1_high <- mean_score_RDA1 + 3 * sd_score_RDA1
#threshold_RDA1_low <- mean_score_RDA1 - 3 * sd_score_RDA1
#top10_outliers_RDA1 <- head(lociscores[order(lociscores$RDA1, decreasing = TRUE), ],10)
#lowest10_outliers_RDA1 <- head(lociscores[order(lociscores$RDA1, decreasing = FALSE), ],10)
#
##sorted_outliers_RDA2 <- lociscores[order(lociscores$RDA1, decreasing = TRUE), ]
#mean_score_RDA2 <- mean(lociscores$RDA2)
#sd_score_RDA2 <- sd(lociscores$RDA2)
#threshold_RDA2 <- mean_score + 3 * sd_score
#threshold_RDA2_high <- mean_score_RDA2 + 3 * sd_score_RDA2
#threshold_RDA2_low <- mean_score_RDA2 - 3 * sd_score_RDA2
#top10_outliers_RDA2 <- head(lociscores[order(lociscores$RDA2, decreasing = TRUE), ],10)
#
#
## Outliers on the posit^ive side of RDA1
#positive_outliers_RDA1 <- top10_outliers_RDA1[which(top10_outliers_RDA1$RDA1 > threshold_RDA1_high),]
#negative_outliers_RDA1 <- lowest10_outliers_RDA1[which(lowest10_outliers_RDA1$RDA1 < threshold_RDA1_high),]
#
#positive_outliers_RDA2 <- top10_outliers_RDA2[which(top10_outliers_RDA2$RDA2 > threshold_RDA2_high),]
## Outliers on the negative side of RDA1
##negative_outliers <- which(top10_outliers_RDA1$RDA1 < threshold_RDA1_low)
## Outliers on the positive side of RDA1
#plot_list <- list()
##positive_outliers <- lociscores$loci[lociscores$RDA1 > threshold ]
#for (i in 1:nrow(positive_outliers_RDA1)) {
#  locus <- positive_outliers_RDA1$loci[i]
#  # Extract the relevant columns from data.frame1 and data.frame2
#  AlleleFrequency <- AllFreq[locus][,1]
#  Bio9 <- Env2$wc2.1_2.5m_bio_9
#  correlations <- data.frame(AlleleFrequency = AlleleFrequency, Bio9 = Bio9)
#  # Calculate correlation coefficient
#  correlation_coef <- cor(AlleleFrequency, Bio9)
#  # Create ggplot object for each locus
#  p <- ggplot(correlations, aes(x = Bio9, y = AlleleFrequency)) +
#    geom_point() +
#    geom_smooth(method = "lm") +
#    labs(title = paste("Locus:", locus, "- Correlation:", round(correlation_coef, 2))) +
#    theme_minimal()
#  # Store ggplot object in the list
#  plot_list[[i]] <- p
#}
#
## Combine all ggplot objects into one plot with multiple grids (facets)
#multi_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
#
#png(filename = paste0(outdir,"/MultiPlot_Correlations.png"), width = 800, height = 600)
#print(multi_plot)
#dev.off()
#
#
#
######NEGATIVE OUTLIERS RDA1
#plot_list <- list()
##positive_outliers <- lociscores$loci[lociscores$RDA1 > threshold ]
#for (i in 1:nrow(negative_outliers_RDA1)) {
#  locus <- negative_outliers_RDA1$loci[i]
#  # Extract the relevant columns from data.frame1 and data.frame2
#  AlleleFrequency <- AllFreq[locus][,1]
#  Bio9 <- Env2$wc2.1_2.5m_bio_9
#  correlations <- data.frame(AlleleFrequency = AlleleFrequency, Bio9 = Bio9)
#  # Calculate correlation coefficient
#  correlation_coef <- cor(AlleleFrequency, Bio9)
#  # Create ggplot object for each locus
#  p <- ggplot(correlations, aes(x = Bio9, y = AlleleFrequency)) +
#    geom_point() +
#    geom_smooth(method = "lm") +
#    labs(title = paste("Locus:", locus, "- Correlation:", round(correlation_coef, 2))) +
#    theme_minimal()
#  # Store ggplot object in the list
#  plot_list[[i]] <- p
#}
#
## Combine all ggplot objects into one plot with multiple grids (facets)
#multi_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
#print(multi_plot)
#
#
#
####now compare with RDA2 to "cross validate"
#plot_list <- list()
##positive_outliers <- lociscores$loci[lociscores$RDA1 > threshold ]
#for (i in 1:nrow(positive_outliers_RDA2)) {
#  locus <- positive_outliers_RDA2$loci[i]
#  # Extract the relevant columns from data.frame1 and data.frame2
#  AlleleFrequency <- AllFreq[locus][,1]
#  Bio9 <- Env2$wc2.1_2.5m_bio_9
#  correlations <- data.frame(AlleleFrequency = AlleleFrequency, Bio9 = Bio9)
#  # Calculate correlation coefficient
#  correlation_coef <- cor(AlleleFrequency, Bio9)
#  # Create ggplot object for each locus
#  p <- ggplot(correlations, aes(x = Bio9, y = AlleleFrequency)) +
#    geom_point() +
#    geom_smooth(method = "lm") +
#    labs(title = paste("Locus:", locus, "- Correlation:", round(correlation_coef, 2))) +
#    theme_minimal()
#  # Store ggplot object in the list
#  plot_list[[i]] <- p
#}
#
## Combine all ggplot objects into one plot with multiple grids (facets)
#multi_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
#print(multi_plot)
#
#
####now compare with other bioclim to "cross validate"
#plot_list <- list()
##positive_outliers <- lociscores$loci[lociscores$RDA1 > threshold ]
#for (i in 1:nrow(positive_outliers_RDA1)) {
#  locus <- positive_outliers_RDA1$loci[i]
#  # Extract the relevant columns from data.frame1 and data.frame2
#  AlleleFrequency <- AllFreq[locus][,1]
#  Bio4 <- Env2$wc2.1_2.5m_bio_4
#  correlations <- data.frame(AlleleFrequency = AlleleFrequency, Bio4 = Bio4)
#  # Calculate correlation coefficient
#  correlation_coef <- cor(AlleleFrequency, Bio4)
#  # Create ggplot object for each locus
#  p <- ggplot(correlations, aes(x = Bio4, y = AlleleFrequency)) +
#    geom_point() +
#    geom_smooth(method = "lm") +
#    labs(title = paste("Locus:", locus, "- Correlation:", round(correlation_coef, 2))) +
#    theme_minimal()
#  # Store ggplot object in the list
#  plot_list[[i]] <- p
#}
#multi_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
#print(multi_plot)
#
####RDA 2 with Bio??
#
#plot_list <- list()
##positive_outliers <- lociscores$loci[lociscores$RDA1 > threshold ]
#for (i in 1:nrow(positive_outliers_RDA2)) {
#  locus <- positive_outliers_RDA2$loci[i]
#  # Extract the relevant columns from data.frame1 and data.frame2
#  AlleleFrequency <- AllFreq[locus][,1]
#  Bio4 <- Env2$wc2.1_2.5m_bio_4
#  correlations <- data.frame(AlleleFrequency = AlleleFrequency, Bio4 = Bio4)
#  # Calculate correlation coefficient
#  correlation_coef <- cor(AlleleFrequency, Bio4)
#  # Create ggplot object for each locus
#  p <- ggplot(correlations, aes(x = Bio4, y = AlleleFrequency)) +
#    geom_point() +
#    geom_smooth(method = "lm") +
#    labs(title = paste("Locus:", locus, "- Correlation:", round(correlation_coef, 2))) +
#    theme_minimal()
#  # Store ggplot object in the list
#  plot_list[[i]] <- p
#}
#multi_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
#print(multi_plot)
#
####now compare with other bioclim to "cross validate"
#plot_list <- list()
##positive_outliers <- lociscores$loci[lociscores$RDA1 > threshold ]
#for (i in 1:nrow(positive_outliers)) {
#  locus <- positive_outliers$loci[i]
#  # Extract the relevant columns from data.frame1 and data.frame2
#  AlleleFrequency <- AllFreq[locus][,1]
#  Bio4 <- Env2$wc2.1_2.5m_bio_4
#  correlations <- data.frame(AlleleFrequency = AlleleFrequency, Bio4 = Bio4)
#  # Calculate correlation coefficient
#  correlation_coef <- cor(AlleleFrequency, Bio4)
#  # Create ggplot object for each locus
#  p <- ggplot(correlations, aes(x = Bio4, y = AlleleFrequency)) +
#    geom_point() +
#    geom_smooth(method = "lm") +
#    labs(title = paste("Locus:", locus, "- Correlation:", round(correlation_coef, 2))) +
#    theme_minimal()
#  # Store ggplot object in the list
#  plot_list[[i]] <- p
#}
#
##multi_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
##print(multi_plot)
#
#
#