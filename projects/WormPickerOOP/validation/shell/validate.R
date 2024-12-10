library(raster)
library(rasterVis)
library(rgdal)
library(tidyverse)
library(mapproj)
library(ggpubr)
library(ncdf4)

mkdir "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/validation/bboxAUT/eu_demography/JRC"
cd "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/validation/bboxAUT/eu_demography/JRC"
unzip "/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/validation/bboxAUT/eu_demography/JRC_GRID_2018.zip"

library(terra)
library(tidyverse)


## demography
Original <- rast("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/validation/bboxAUT/eu_demography/JRC/JRC_1K_POP_2018.tif")
Rasdaman <- rast("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/validation/bboxAUT/eu_demography/eu_demography_rasdaman.tiff")

Demography <- read.csv("/media/inter/ssteindl/FC/usecaserepo/SYNC0524/uc3-drosophola-genetics/projects/WormPickerOOP/validation/bboxAUT/eu_demography/eu_dem.csv", 
        na.string = c("-1.5", "-1", "-2", "NA"))
# Load the CSV file with point data

# Convert the data frame to a SpatVector (spatial vector object) with latitude and longitude
points_vect <- vect(Demography, geom = c("lon", "lat"), crs = crs(Original.wgs84))
extracted_values <- terra::extract(Original,points_vect)
combined_data <- cbind(Demography, extracted_values)
cor_results <- cor(Pest[, names(raster_stack)], extracted_values[, -1], use = "pairwise.complete.obs")

# Convert the correlation matrix to a data frame for plotting
cor_df <- na.omit(as.data.frame(as.table(cor_results)))
colnames(cor_df) <- c("Band", "Metric", "Correlation")

# Filter out self-correlations (diagonal) if needed
cor_df <- cor_df[cor_df$Band == cor_df$Metric, ]


# Dot plot of the correlations with vertical x-axis labels
ggplot(cor_df, aes(x = Band, y = Correlation)) +
    geom_point(size = 5, color = "blue") +
    geom_text(aes(label = round(Correlation, 2)), vjust = -1) +
    theme_minimal() +
    ylim(-1, 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
        title = "Correlation between Raster Bands and Point Data",
        x = "Raster Band",
        y = "Correlation Coefficient"
    )

ggsave("/media/inter/mkapun/projects/uc3-drosophola-genetics_old/projects/PlotRaster/RasterMaps_Pest_corr.png",
    width = 12,
    height = 8
)
