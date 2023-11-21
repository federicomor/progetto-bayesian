library(maps)
library(ggplot2)
library(sf)
library(lubridate)
#devtools::install_github("dgrtwo/gganimate")
library(gganimate)
library(av)
library(ggpubr)
library(jpeg)

load("../data/data_agc_lomb_part1.RData")
load("../data/data_agc_lomb_part2.RData")
AGC_Dataset = rbind(parte1, parte2)

head(AGC_Dataset)
data_agc_lomb=AGC_Dataset
data_agc_lomb$Time = as.Date(data_agc_lomb$Time)
df_agri$Time = as.Date(df_agri$Time)

rm(parte1) # save memory
rm(parte2) # save memory
rm(AGC_Dataset) # save memory

sites <- data.frame(
	longitude = unique(df_agri$Longitude), 
	latitude = unique(df_agri$Latitude))

stations = unique(df_agri$IDStations)
cols = colora(6,970,show=F)
pad <- 0.2 * c(-1, 1)
DATE_FORMAT = "%Y-%m-%d"


# carica file shp
details <- 3
details_altre_regioni <- 2

confini  <- c("Emilia-Romagna","Piemonte","Lombardia","Trentino-Alto Adige","Veneto")
regioni_italiane <- st_read(paste0("italia/gadm40_ITA_",details,".shp"))
regioni_italiane_2 <- st_read(paste0("italia/gadm40_ITA_",details_altre_regioni,".shp"))

lombardia <- regioni_italiane[regioni_italiane$NAME_1 == "Lombardia",]
altre_regioni <- regioni_italiane_2[regioni_italiane_2$NAME_1 %in% confini,]

shp_map <- st_read(paste0("italia/gadm40_ITA_",1,".shp"))

