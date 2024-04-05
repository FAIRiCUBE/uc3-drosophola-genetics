library(raster)
library(rasterVis)
library(rgdal)
library(tidyverse)
library(mapproj)
library(ggpubr)
library(ncdf4)

world <- map_data("world")
biod <- getData("worldclim", var = "bio", res = 2.5)
bio1 <- as.data.frame(biod$bio1 / 10, xy = TRUE, na.rm = TRUE)
rownames(bio1) <- c()

setwd("/media/inter/mkapun/projects/uc3-drosophola-genetics/projects/PlotRaster")

META <- read.table("dest_v2.samps_8Jun2023.csv", header = T, sep = ",")
Europe.meta <- META %>%
    filter(continent == "Europe" | continent == "Asia") %>%
    mutate(
        Lat = round(lat, 1),
        Long = round(long, 1)
    ) %>%
    select(Lat, Long) %>%
    group_by(Lat, Long) %>%
    summarize(Count = n())

NorthAmerica.meta <- META %>%
    filter(continent == "North_America") %>%
    mutate(
        Lat = round(lat, 1),
        Long = round(long, 1)
    ) %>%
    select(Lat, Long) %>%
    group_by(Lat, Long) %>%
    summarize(Count = n())


Europe <- ggplot() +
    geom_tile(
        data = bio1,
        aes(fill = bio1, x = x, y = y),
        alpha = 1,
    ) +
    geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        colour = "black",
        fill = NA
    ) +
    geom_point(
        data = Europe.meta,
        aes(
            x = Long,
            y = Lat,
            size = Count
        )
    ) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_fixed(xlim = c(-30, 60), ylim = c(30, 80)) +
    scale_fill_gradientn(
        name = "Temp (°C)",
        colours = c("#0094D1", "#7ec9a5", "#FEED99", "#AF3301"),
        breaks = c(-20, 0, 20, 40)
    ) +
    ggtitle("Europe - Average annual temperature")

NorthAmerica <- ggplot() +
    geom_tile(
        data = bio1,
        aes(fill = bio1, x = x, y = y),
        alpha = 1,
    ) +
    geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        colour = "black",
        fill = NA
    ) +
    geom_point(
        data = NorthAmerica.meta,
        aes(
            x = Long,
            y = Lat,
            size = Count
        )
    ) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_fixed(xlim = c(-165, -40), ylim = c(0, 80)) +
    scale_fill_gradientn(
        name = "Temp (°C)",
        colours = c("#0094D1", "#7ec9a5", "#FEED99", "#AF3301"),
        breaks = c(-20, 0, 20, 40)
    ) +
    ggtitle("North America - Average annual temperature")

PLOT <- ggarrange(Europe, NorthAmerica,
    common.legend = TRUE,
    heights = c(4, 4),
    ncol = 2,
    align = "h",
    legend = "bottom"
)
# ggsave("Map.pdf",
#     PLOT,
#     width = 12,
#     height = 5
# )
ggsave("Map.png",
    PLOT,
    width = 12,
    height = 5
)


### Do the same with netcdf data

nc_data <- nc_open("NSAT20150909_timeless.nc")
# Save the print(nc) dump to a text file
{
    sink("NSAT20150909_timeless.txt")
    print(nc_data)
    sink()
}

lon <- -1 * ncvar_get(nc_data, "lon")
lat <- -1 * ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")

ndvi.array <- ncvar_get(nc_data, "Tair") # store the data in a 3-dimensional array
dim(ndvi.array)

fillvalue <- ncatt_get(nc_data, "Tair", "_FillValue")
fillvalue

ndvi.array[ndvi.array == fillvalue$value] <- NA

r <- raster(t(ndvi.array), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
bio_NSAT <- as.data.frame(r$layer / 10, xy = TRUE, na.rm = TRUE)

Europe_bio <- ggplot() +
    geom_tile(
        data = bio_NSAT,
        aes(fill = layer, x = x, y = -y),
        alpha = 1,
    ) +
    geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        colour = "black",
        fill = NA
    ) +
    geom_point(
        data = Europe.meta,
        aes(
            x = Long,
            y = Lat,
            size = Count
        )
    ) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_fixed(xlim = c(-30, 60), ylim = c(30, 80)) +
    scale_fill_gradientn(
        name = "Temp (°C)",
        colours = c("#0094D1", "#7ec9a5", "#FEED99", "#AF3301"),
        breaks = c(-20, 0, 20, 40)
    ) +
    ggtitle("Europe")

NorthAmerica_bio <- ggplot() +
    geom_tile(
        data = bio_NSAT,
        aes(fill = layer, x = x, y = -y),
        alpha = 1,
    ) +
    geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        colour = "black",
        fill = NA
    ) +
    geom_point(
        data = NorthAmerica.meta,
        aes(
            x = Long,
            y = Lat,
            size = Count
        )
    ) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    coord_fixed(xlim = c(-165, -40), ylim = c(0, 80)) +
    scale_fill_gradientn(
        name = "Temp (°C)",
        colours = c("#0094D1", "#7ec9a5", "#FEED99", "#AF3301"),
        breaks = c(-20, 0, 10, 40)
    ) +
    ggtitle("North America - ")
NorthAmerica_bio

PLOT <- ggarrange(Europe, NorthAmerica, Europe_bio, NorthAmerica_bio,
    common.legend = TRUE,
    heights = c(4, 4),
    ncol = 2,
    nrow = 2,
    align = "h",
    legend = "bottom"
)
# ggsave("Map_NSAT.pdf",
#     PLOT,
#     width = 12,
#     height = 10
# )
ggsave("Map_NSAT.png",
    PLOT,
    width = 12,
    height = 10
)

##
bio <- raster::extract(brick("NSAT20150909_timeless.nc"), Europe.meta[, c(2, 1)])
Europe.bio <- cbind(Europe.meta, bio[, 1] / 10)
colnames(Europe.bio) <- c("Lat", "Long", "Count", "Temp")

Europe_plot <- ggplot(Europe.bio, aes(x = Long, y = Lat, fill = Temp)) +
    geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        colour = "black",
        fill = "grey95"
    ) +
    geom_point(
        aes(
            size = Count
        ),
        pch = 21, colour = "black"
    ) +
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    scale_fill_gradientn(
        name = "Temp (°C)",
        colours = c("#0094D1", "#7ec9a5", "#FEED99", "#AF3301"),
    ) +
    coord_fixed(xlim = c(-30, 60), ylim = c(30, 80)) +
    ggtitle("Europe")
Europe_plot

ggsave("Europe.png",
    Europe_plot,
    width = 5.5,
    height = 4
)
