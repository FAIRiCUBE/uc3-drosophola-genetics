args <- commandArgs(TRUE)

library(dplyr)
library(ggplot2)
library(grid)

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

meta <- args[1]
outclim<-args[2]

meta_object <- get_meta_data(meta)
worldclimdata <- get_worldclim_data(meta_object)
saveRDS(worldclimdata, file=outclim)
