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
  meta.sub <- meta %>% dplyr::select(sampleId, lat, long)
  #meta.sub$sampleId <- gsub("\\-", ".", meta.sub$sampleId)
  return(meta.sub)
}

get_worldclim_data <- function(meta.sub) {
  biod <- worldclim_global(var = "bio", 2.5, "data")
  bio <- raster::extract(biod, cbind(meta.sub$long, meta.sub$lat))
  bio.sub <- as.data.frame(bio)
  return(bio.sub)
}

meta2024 <- "/home/ssteindl/mounts/BioMem_2/ssteindl/UC3/ClimateData/samples_europe_pass.csv"
meta <- get_meta_data(meta2024)
WC <- get_worldclim_data(meta)
rownames(WC) <- meta$sampleId

Env <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/ClimateData/Env.csv", header = TRUE)
