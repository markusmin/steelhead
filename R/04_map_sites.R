##### 04_map_sites #####

# Description: This will map the detection sites 

library(tidyverse)
# library(readxl)
library(here)
library(rgdal)
library(broom)
library(ggrepel)
# library(ggpubr)
# library(sf)
# library(lubridate)


##### Load shapefiles ##### 

# USA boundaries
usa_spdf <- readOGR(dsn = here::here("map_files", "USA_adm0.shp"))
usa_spdf_fort <- tidy(usa_spdf)
# 
# # Major rivers
rivers_spdf <- readOGR(dsn = here::here("map_files", "NA_Lakes_and_Rivers","hydrography_l_rivers_v2.shp"))
summary(rivers_spdf)
# proj4string(rivers_spdf)
# # rivers_spdf_transform <- spTransform(rivers_spdf, CRS(usa_spdf))
rivers_spdf_transform <- spTransform(rivers_spdf, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
# # rivers_spdf_transform <- coordinates(spTransform(rivers_spdf, CRS("+proj=laea +datum=WGS84")))
# # rivers_spdf_fort <- tidy(rivers_spdf)
rivers_spdf_fort <- tidy(rivers_spdf_transform)

# CRB streams
# Data from here: https://www.fisheries.noaa.gov/resource/data/columbia-basin-historical-ecology-project-data
CRB_streams_spdf <- readOGR(dsn = here::here("map_files", "crb_streams_over8m_100k","crb_streams_over8m_100k.shp"))
CRB_streams_spdf_transform <- spTransform(CRB_streams_spdf, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# We have to subset this to the rivers we want
CRB_streams_spdf_transform$y_coord
CRB_stream_names <-  data.frame(unique(CRB_streams_spdf_transform$GNIS_NAME))
colnames(CRB_stream_names) <- c("GNIS_NAME")
# Streams we need: Columbia River, Hood River, Fifteenmile Creek, Deschutes River, 
# John Day River, Umatilla River, Yakima River, Wenatchee River, Entiat River,
# Walla Walla River, Tucannon River, Clearwater River, Grande Ronde River,
# Imnaha River, Snake River, Salmon River

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
                    CRB_stream_names$GNIS_NAME[grep("Salmon", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Bridge", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("White Creek", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Catherine", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Looking", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Lapwai", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Klickitat", CRB_stream_names$GNIS_NAME)],
                    CRB_stream_names$GNIS_NAME[grep("Trout", CRB_stream_names$GNIS_NAME)])


# Remove the incorrectly selected ones
incorrect_streams <- c("Hoodoo Creek", "Little Wenatchee River", 
                       "Clearwater Creek", "Little Clearwater River", "Snake Creek", 
                       "No Snake Creek", "White Salmon River", "Salmon Falls Creek",
                       "Salmon Creek", "Little White Salmon River", "Little Salmon River",
                       "South Fork Salmon Falls Creek", "Little Salmon Creek",
                       "Big Salmon Creek", "Salmon la Sac Creek", "Trout Lake Creek")

subset_streams <- subset_streams[!(subset_streams %in% incorrect_streams)]

CRB_streams_spdf_transform_subset <- subset(CRB_streams_spdf_transform, GNIS_NAME %in% subset_streams)

# proj4string(CRB_streams_spdf)
CRB_streams_fort <- tidy(CRB_streams_spdf_transform_subset)
# Remove shapefiles to save space in workspace
# rm(CRB_streams_spdf_transform_subset, CRB_streams_spdf_transform)

##### Intermediate maps/subsetting #####

## Figure out which of these rivers we need to subset
# NOTE: Need to clip out the coastline.

subset(rivers_spdf_fort, lat < 46.2 & lat > 46.15 & long < -123.10 & long > -123.35 |
         lat < 46 & lat > 45.95 & long < -122.72 & long > -122.9 |
         lat < 45.61 & lat > 45.5 & long < -122.16 & long > -122.31 |
         lat < 45.75 & lat > 45.68 & long < -120.67 & long > -120.9 |
         lat < 45.8 & lat > 45.6 & long < -120.1 & long > -120.2 |
         lat < 46.1 & lat > 46 & long < -118.9 & long > -118.95 |
         lat < 46 & lat > 45.9 & long < -119.2 & long > -119.3 |
         lat < 46.48 & lat > 46.35 & long < -119.2 & long > -119.3 |
         lat < 46.3 & lat > 46.15 & long < -118.9 & long > -119.1 |
         lat < 47 & lat > 46.9 & long < -119.95 & long > -120.15 |
         lat < 48.08 & lat > 48 & long < -119.6 & long > -120 |
         lat < 47.95 & lat > 47.9 & long < -118.5 & long > -118.8 |
         lat < 48.4 & lat > 48.2 & long < -117.9 & long > -118.4 |
         lat < 46.75 & lat > 46.6 & long < -117.35 & long > -117.8 |
         lat < 46.67 & lat > 46.6 & long < -116 & long > -116.33) -> columbia_fragments
subset(rivers_spdf_fort, id %in% unique(columbia_fragments$id)) -> rivers_subset

CRB_map <- ggplot(usa_spdf_fort, aes(x = long, y = lat, group = group))+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_fixed(ylim = c(44.3,48.45),  xlim = c(-125.15,-113.23), ratio = 1.3)+
  # Polygons for base map
  # geom_polygon(color = "gray70", size = 0.2, fill = rgb(251, 234, 194, max=255))+
  # Polygons for CRB streams
  # geom_polygon(data = CRB_streams_fort, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "gray70", size = 0.2, fill = "gray70") +
  # Lines for North American Rivers (includes Columbia and Snake)
  geom_path(data = rivers_subset, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "gray70", size = 0.2, fill = "gray70") +
  # geom_path(data = CRB_streams_fort, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "black", size = 0.3) +
  # Polygon for mainstem
  # geom_polygon(data = CRB_mainstem_spdf_transform_fort, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "gray70", size = 0.2, fill = "gray70") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  guides(fill = guide_legend(title = "Legend"))

# ggsave(here("figures", "CRB_subset_map.pdf"), CRB_map, height = 6, width  = 10)  



##### MAKE BASE MAP #####
CRB_map <- ggplot(usa_spdf_fort, aes(x = long, y = lat, group = group))+
  ylab("Latitude")+
  xlab("Longitude")+
  coord_fixed(ylim = c(44.3,48.45),  xlim = c(-124.5,-114), ratio = 1.3)+
  # Polygons for base map
  geom_polygon(color = "gray70", size = 0.2, fill = rgb(251, 234, 194, max=255))+
  # Polygons for CRB streams
  geom_polygon(data = CRB_streams_fort, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "gray70", size = 0.2, fill = "gray70") +
  # Lines for North American Rivers (includes Columbia and Snake)
  geom_path(data = rivers_subset, aes(x = long, y = lat, group = group), inherit.aes = FALSE, color = "gray70", size = 0.5, fill = "gray70") +
  # ADD LABELS FOR RIVERS
  # annotate("segment", x = -122.84, y = 47.68, xend = -123.18, yend = 47.68, size = 0.5, lty = 1) + # Columbia River
  annotate("text", x = -123.33, y = 46.27, label = "Columbia River", size = 3.5, fontface = 'italic', hjust = 0) + # Columbia River
  annotate("text", x = -121.7, y = 45.5, label = "Hood R.", size = 2.5, fontface = 'italic', hjust = 1) + # Hood River
  annotate("text", x = -121.2, y = 45.44, label = "Fifteenmile Cr.", size = 2.5, fontface = 'italic', hjust = 1, angle = 45) + # Fifteenmile Cr.
  annotate("text", x = -120.94, y = 44.8, label = "Deschutes R.", size = 2.5, fontface = 'italic', hjust = 1) + # Deschutes R.
  annotate("text", x = -119.8, y = 44.92, label = "John Day R.", size = 2.5, fontface = 'italic', hjust = 1) + # John Day R.
  annotate("text", x = -118.98, y = 45.60, label = "Umatilla R.", size = 2.5, fontface = 'italic', hjust = 1) + # Umatilla R.
  annotate("text", x = -118, y = 46.15, label = "Walla Walla R.", size = 2.5, fontface = 'italic', hjust = 1) + # Walla Walla R.
  annotate("text", x = -117.94, y = 46.396, label = "Tucannon R.", size = 2.5, fontface = 'italic', hjust = 1) + # Tucannon R.
  annotate("text", x = -117.3, y = 45.70, label = "Grande Ronde R.", size = 2.5, fontface = 'italic', hjust = 1) + # Grande Ronde R.
  annotate("text", x = -115.4, y = 46.475, label = "Clearwater R.", size = 2.5, fontface = 'italic', hjust = 1) + # Clearwater R.
  annotate("text", x = -116.95, y = 45.52, label = "Imanha R.", size = 2.5, fontface = 'italic', hjust = 1) + # Imanha R.
  annotate("text", x = -116.579, y = 44.997, label = "Snake R.", size = 3, fontface = 'italic', hjust = 1) + # Snake R.
  annotate("text", x = -115.05, y = 45.65, label = "Salmon R.", size = 2.5, fontface = 'italic', hjust = 1) + # Salmon R.
  annotate("text", x = -120.44, y = 46.40, label = "Yakima R.", size = 2.5, fontface = 'italic', hjust = 1) + # Yakima R.
  annotate("text", x = -120.62, y = 47.50, label = "Wenatchee R.", size = 2.5, fontface = 'italic', hjust = 1) + # Wenatchee R.
  annotate("text", x = -120.33, y = 47.88, label = "Entiat R.", size = 2.5, fontface = 'italic', hjust = 1) + # Entiat R.
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

##### Plot JDR sites on top of base map #####
JDR_event_det_counts <- read.csv(here::here("model_files", "JDR_event_det_counts.csv"))


JDR_site_map <- CRB_map +
  # Add dams
  geom_point(data = subset(JDR_event_det_counts, dam == "dam"), 
             aes(x = event_site_longitude, y = event_site_latitude), 
             shape = 73, size = 4, inherit.aes = FALSE) +
  geom_text(data = subset(JDR_event_det_counts, dam == "dam"), 
            aes(x = event_site_longitude, y = event_site_latitude+0.2, label = dam_abbr), 
            size = 3, inherit.aes = FALSE) +
  # Add detection sites
  geom_point(data = JDR_event_det_counts, aes(x = event_site_longitude, 
                                              y = event_site_latitude), size = 1, inherit.aes = FALSE)
  # geom_text_repel(data = JDR_event_det_counts, aes(x = event_site_longitude,
  #                                                  y = event_site_latitude+0.02, label = event_site_name), 
  #                 max.overlaps = 100, size = 1, inherit.aes = FALSE)

ggsave(here("figures", "JDR_site_map2.pdf"), JDR_site_map, height = 6, width  = 10)  