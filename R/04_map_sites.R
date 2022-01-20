##### 04_map_sites #####

# Description: This will map the detection sites 

library(tidyverse)
# library(readxl)
library(here)
library(rgdal)
library(broom)
# library(ggpubr)
# library(sf)
# library(lubridate)


##### Load shapefiles ##### 

# USA boundaries
usa_spdf <- readOGR(dsn = here("map_files", "USA_adm0.shp"))
usa_spdf_fort <- tidy(usa_spdf)
# 
# # Major rivers
# rivers_spdf <- readOGR(dsn = here("map_files", "NA_Lakes_and_Rivers","hydrography_l_rivers_v2.shp"))
# summary(rivers_spdf)
# proj4string(rivers_spdf)
# # rivers_spdf_transform <- spTransform(rivers_spdf, CRS(usa_spdf))
# rivers_spdf_transform <- spTransform(rivers_spdf, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
# # rivers_spdf_transform <- coordinates(spTransform(rivers_spdf, CRS("+proj=laea +datum=WGS84")))
# # rivers_spdf_fort <- tidy(rivers_spdf)
# rivers_spdf_fort <- tidy(rivers_spdf_transform)

# CRB streams
# Data from here: https://www.fisheries.noaa.gov/resource/data/columbia-basin-historical-ecology-project-data
CRB_streams_spdf <- readOGR(dsn = here("map_files", "crb_streams_over8m_100k","crb_streams_over8m_100k.shp"))
CRB_streams_spdf_transform <- spTransform(CRB_streams_spdf, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# We have to subset this to the rivers we want
CRB_streams_spdf_transform$y_coord
CRB_stream_names <-  data.frame(unique(CRB_streams_spdf_transform$GNIS_NAME))
colnames(CRB_stream_names) <- c("GNIS_NAME")
# Streams we need: Columbia River, Hood River, Fifteenmile Creek, Deschutes River, 
# John Day River, Umatilla River, Yakima River, Wenatchee River, Entiat River,
# Walla Walla River, Tucannon River, Clearwater River, Grande Ronde River,
# Imnaha River, Snake River, Salmon River

CRB_stream_names$GNIS_NAME[grep("", CRB_stream_names$GNIS_NAME)]

subset_streams <- c(CRB_stream_names$GNIS_NAME[grep("Columbia", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Hood", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Fifteenmile", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Deschutes", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("John Day", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Umatilla", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Yakima", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Wenatchee", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Entiat", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Walla Walla", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Tucannon", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Clearwater", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Grande Ronde", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Imnaha", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Snake", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Salmon", CRB_stream_names$GNIS_NAME)])

# Remove the incorrectly selected ones
incorrect_streams <- c("Hoodoo Creek", "Little Deschutes River", "Little Wenatchee River", 
                       "Clearwater Creek", "Little Clearwater River", "Snake Creek", 
                       "No Snake Creek", "White Salmon River", "Salmon Falls Creek",
                       "Salmon Creek", "Little White Salmon River", "Little Salmon River",
                       "South Fork Salmon Falls Creek", "Little Salmon Creek",
                       "Big Salmon Creek", "Salmon la Sac Creek")

subset_streams <- subset_streams[!(subset_streams %in% incorrect_streams)]

CRB_streams_spdf_transform_subset <- subset(CRB_streams_spdf_transform, GNIS_NAME %in% subset_streams)

# proj4string(CRB_streams_spdf)
CRB_streams_fort <- tidy(CRB_streams_spdf_transform_subset)

##### Mainstem CRB
CRB_mainstem_spdf <- readOGR(dsn = here("map_files", "crb_habitat_model_over8m_100k","crb_habitat_model_over8m_100k.shp"))
CRB_mainstem_spdf_transform <- spTransform(CRB_mainstem_spdf, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
CRB_mainstem_spdf_transform_fort <- tidy(CRB_mainstem_spdf_transform)
CRB_mainstem_spdf_transform$
##### Make base map #####

CRB_map <- ggplot(usa_spdf_fort, aes(x = long, y = lat, group = group))+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_fixed(ylim = c(44.3,48.45),  xlim = c(-125.15,-113.23), ratio = 1.3)+
  # Polygons for base map
  geom_polygon(color = "gray70", size = 0.2, fill = rgb(251, 234, 194, max=255))+
  # Polygons for river
  geom_polygon(data = CRB_streams_fort, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "gray70", size = 0.2, fill = "gray70") +
  # geom_path(data = CRB_streams_fort, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "black", size = 0.3) +
  # Polygon for mainstem
  geom_polygon(data = CRB_mainstem_spdf_transform_fort, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "gray70", size = 0.2, fill = "gray70") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  guides(fill = guide_legend(title = "Legend"))
# scale_x_continuous(breaks = seq(-128, -123, 1), expand = c(0,0), labels=c(expression(paste(128*degree,"W")),
#                                                                           expression(paste(127*degree,"W")),
#                                                                           expression(paste(126*degree,"W")),
#                                                                           expression(paste(125*degree,"W")),
#                                                                           expression(paste(124*degree,"W")),
#                                                                           expression(paste(123*degree,"W"))))+
# scale_y_continuous(breaks = seq(47, 51, 1), expand = c(0,0), labels=c(expression(paste(47*degree,"N")),
#                                                                       expression(paste(48*degree,"N")),
#                                                                       expression(paste(49*degree,"N")),
#                                                                       expression(paste(50*degree,"N")),
#                                                                       expression(paste(51*degree,"N"))))+


ggsave(here("figures", "CRB_map_v3.pdf"), CRB_map, height = 6, width  = 10)  

