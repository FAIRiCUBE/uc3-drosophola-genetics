library(raster)
# library(rasterVis)
# library(rgdal)
library(tidyverse)
# library(mapproj)
# library(ggpubr)
library(ncdf4)
library(ggmap)
library(osmdata)
# install.packages("osmdata")
# install.packages("rgdal", repos="http://R-Forge.R-project.org")
# install.packages("OpenStreetMap")
library(sp)

### get data from Vienna
Vienna <- get_stadiamap(getbb("Vienna"), source = "stadia", maptype = "stamen_toner", zoom = 13)

### read WorldClim data
world <- map_data("world")
biod <- getData("worldclim", var = "bio", res = 2.5)
bio1 <- as.data.frame(biod$bio1 / 10, xy = TRUE, na.rm = TRUE)
rownames(bio1) <- c()

### read data with participants
DATA <- read.csv("/media/inter/mkapun/projects/uc3-drosophola-genetics/projects/ViennaCityFly2004/VCF2024_Participants_Geocode.csv", header = T)

### Read NetCDF data from GeoSphere
nc_data <- nc_open("/media/inter/mkapun/projects/uc3-drosophola-genetics/projects/ViennaCityFly2004/SPARTACUS - Spatial Dataset for Climate in Austria Datensatz_202303_202303.nc")

## dump metainfo to text file
{
    sink("/media/inter/mkapun/projects/uc3-drosophola-genetics/projects/ViennaCityFly2004/SPARTACUS - Spatial Dataset for Climate in Austria Datensatz_202303_202303.txt")
    print(nc_data)
    sink()
}

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")

ndvi.array <- ncvar_get(nc_data, "TM") # store the data in a 3-dimensional array
dim(ndvi.array)

fillvalue <- ncatt_get(nc_data, "TM", "_FillValue")
fillvalue

ndvi.array[ndvi.array == fillvalue$value] <- NA

r <- raster(t(ndvi.array), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

midatlas <- extent(c(min(lon), max(lon), min(lat), max(lat)))
ma.f <- raster(ext = midatlas, ncol = 500, nrow = 500)

ma.r <- resample(x = r, y = ma.f, method = "bilinear")

bio_NSAT <- as.data.frame(ma.r$layer / 10, xy = TRUE, na.rm = TRUE)

MAP <- ggmap(Vienna) +
    geom_tile(
        data = bio_NSAT,
        aes(fill = layer, x = x, y = y),
        alpha = 0.5
    ) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    scale_fill_gradientn(
        name = "Temp (°C)",
        colours = c("#0094D1", "#7ec9a5", "#FEED99", "#AF3301"),
        breaks = c(-20, 0, 20, 40)
    ) +
    geom_point(
        data = DATA,
        aes(
            x = lon,
            y = lat,
            size = 1
        )
    ) +
    ggtitle("ViennaCityFly 2024") + theme(legend.position = "none")

ggsave(
    file = "/media/inter/mkapun/projects/uc3-drosophola-genetics/projects/ViennaCityFly2004/MAP.pdf",
    MAP,
    width = 8,
    height = 6
)

ggsave(
    file = "/media/inter/mkapun/projects/uc3-drosophola-genetics/projects/ViennaCityFly2004/MAP.png",
    MAP,
    width = 8,
    height = 6
)
