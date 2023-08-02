library(tidyverse)
library(eurostat)
install.packages()
library(leaflet)
library(sf)
library(scales)
library(cowplot)
library(ggthemes)

get_eurostat_geospatial(resolution = 10, 
                        nuts_level = 0, 
                        year = 2016)
