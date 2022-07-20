### 03 - Complete Detection histories - VERSION 2 - FOR HYAK

# For testing:
# setwd("/Users/markusmin/Documents/CBR/steelhead/to_hyak_transfer/2022-07-18-complete_det_hist/")


# This R script adds "implicit site usage" for individuals, i.e., adds in sites
# that must be visited when going from one state to another

# NOTE: To run this script, you will first need to run 02_detection_histories.R
# script, and use the file located here as input: here("model_files", "complete_det_hist.csv")


### Load libraries
# library(here)
library(tidyverse)
library(lubridate)
# library(janitor)

# Read in metadata
complete_event_site_metadata <- read.csv("complete_event_site_metadata.csv", row.names = 1)

##### Assign detection sites to parts of river ####

##### Dams in order - separate Columbia and Snake #####
columbia_dams <- c("Bonneville Adult Fishways (combined)", "McNary Adult Fishways (combined)", 
                   "Priest Rapids Adult Fishways (combined)","Rock Island Adult Fishways (combined)", 
                   "RRF - Rocky Reach Fishway", "Wells Dam Adult Fishways (combined)")
snake_dams <- c("Ice Harbor Adult Fishways (combined)",  "Lower Granite Dam Adult Fishways (combined)")

##### Tributary detection sites

# Natal tributary detection sites
# Here, we need lists of tributary detection sites for each of the natal origins

# Asotin Creek
subset(complete_event_site_metadata, event_site_subbasin_name == "Lower Snake-Asotin") %>%  # asotin creek (Tenmile Creek is not a tributary to Asotin Creek)
  subset(., event_site_name != "TENMC2 - Tenmile Creek, tributary to Snake River") -> ASO_sites_df

# Clearwater River
subset(complete_event_site_metadata, event_site_basin_name == "Clearwater") -> CLE_sites_df

# Deschutes River
subset(complete_event_site_metadata, event_site_basin_name == "Deschutes") -> DES_sites_df

# Entiat River
subset(complete_event_site_metadata, event_site_subbasin_name == "Upper Columbia-Entiat") %>%
  filter(!grepl("Dam", event_site_name)) %>%   # remove all of the dam observations and all of the mainstem observations
  filter(!grepl("RRF", event_site_name)) %>% 
  filter(!grepl("RIA", event_site_name)) %>% 
  filter(!grepl("Rock Island", event_site_name)) %>% 
  filter(!grepl("RIS", event_site_name)) %>% 
  filter(!grepl("EBO", event_site_name)) -> ENT_sites_df

# Fifteenmile Creek
subset(complete_event_site_metadata, event_site_name %in% complete_event_site_metadata$event_site_name[grep("Fifteenmile", complete_event_site_metadata$event_site_name)]) -> FIF_sites_df

# Grande Ronde River
subset(complete_event_site_metadata, event_site_subbasin_name %in% c("Lower Grande Ronde", "Upper Grande Ronde", "Wallowa")) -> GRRO_sites_df

# Hood River
subset(complete_event_site_metadata, event_site_name %in%
         c(complete_event_site_metadata$event_site_name[grep("Hood", complete_event_site_metadata$event_site_name)],
           "EFD - East Fork Diversion Fishway", "PARK - Parkdale Hatchery",
           "MVF - Moving Falls Fish Ladder", "SND - Sandtrap Acclimation Site")) -> HOOD_sites_df

# Imnaha River
subset(complete_event_site_metadata, event_site_subbasin_name == "Imnaha") -> IMN_sites_df

# John Day River
subset(complete_event_site_metadata, event_site_subbasin_name %in% 
         c("Upper John Day", "North Fork John Day", "Middle Fork John Day", "Lower John Day")) -> JDR_sites_df

# Methow River
subset(complete_event_site_metadata, event_site_subbasin_name == "Methow") -> MET_sites_df

# Okanogan River
subset(complete_event_site_metadata, event_site_subbasin_name == "Okanogan") -> OKA_sites_df

# Salmon River
subset(complete_event_site_metadata, event_site_basin_name == "Salmon") -> SAL_sites_df 

# Tucannon River
subset(complete_event_site_metadata, event_site_name %in%
         c(complete_event_site_metadata$event_site_name[grep("Tucannon", complete_event_site_metadata$event_site_name)],
           "CURP (Curl Lake Rearing Pond)")) -> TUC_sites_df

# Umatilla River
subset(complete_event_site_metadata, event_site_subbasin_name == "Umatilla") -> UMA_sites_df

# Walla Walla River
subset(complete_event_site_metadata, event_site_subbasin_name == "Walla Walla") -> WAWA_sites_df

# Wenatchee River
subset(complete_event_site_metadata, event_site_subbasin_name == "Wenatchee") -> WEN_sites_df

subset(complete_event_site_metadata, event_site_basin_name == "Yakima") -> YAK_sites_df # Yakima River

# combine these into one DF to store all of the tributaries that we are interested in - everything else will just be "other tributaries", sorted by reach
ASO_sites_df %>% 
  bind_rows(., CLE_sites_df) %>% 
  bind_rows(., DES_sites_df) %>% 
  bind_rows(., ENT_sites_df) %>% 
  bind_rows(., FIF_sites_df) %>% 
  bind_rows(., GRRO_sites_df) %>% 
  bind_rows(., HOOD_sites_df) %>% 
  bind_rows(., IMN_sites_df) %>% 
  bind_rows(., JDR_sites_df) %>% 
  bind_rows(., MET_sites_df) %>% 
  bind_rows(., OKA_sites_df) %>% 
  bind_rows(., SAL_sites_df) %>% 
  bind_rows(., TUC_sites_df) %>% 
  bind_rows(., UMA_sites_df) %>% 
  bind_rows(., WAWA_sites_df) %>% 
  bind_rows(., WEN_sites_df) %>% 
  bind_rows(., YAK_sites_df) -> origin_sites_df

# get just the site names
# all together
origin_sites <- origin_sites_df$event_site_name

# each individually
ASO_sites <- ASO_sites_df$event_site_name
CLE_sites <- CLE_sites_df$event_site_name
DES_sites <- DES_sites_df$event_site_name
ENT_sites <- ENT_sites_df$event_site_name
FIF_sites <- FIF_sites_df$event_site_name
GRRO_sites <- GRRO_sites_df$event_site_name
HOOD_sites <- HOOD_sites_df$event_site_name
IMN_sites <- IMN_sites_df$event_site_name
JDR_sites <- JDR_sites_df$event_site_name
MET_sites <- MET_sites_df$event_site_name
OKA_sites <- OKA_sites_df$event_site_name
SAL_sites <- SAL_sites_df$event_site_name
TUC_sites <- TUC_sites_df$event_site_name
UMA_sites <- UMA_sites_df$event_site_name
WAWA_sites <- WAWA_sites_df$event_site_name
WEN_sites <- WEN_sites_df$event_site_name
YAK_sites <- YAK_sites_df$event_site_name

# Now we need to figure out which are the sites that are tributaries that we don't care about
# Approach: Subset by basin/subbasin, take out dams and the inriver sites (when applicable), take out the origin sites
# Use the PTAGIS map to figure this out

# Straying sites between BON and MCN:
# subbasins: Middle Columbia-Hood, Umatilla, Middle Columbia-Lake Wallula, Willow, Lower John Day, Lower Deschutes, Trout, Upper Deschutes, 
# Lower Crooked, Upper Crooked, Little Deschutes, Beaver-South Fork, North Fork John Day, Middle Fork John Day, Upper John Day, Klickitat

BON_MCN_other_trib_sites_df <- subset(complete_event_site_metadata, 
                                      event_site_subbasin_name %in% c("Middle Columbia-Hood", "Umatilla", "Middle Columbia-Lake Wallula", "Willow", 
                                                                      "Lower John Day", "Lower Deschutes", "Trout", "Upper Deschutes", 
                                                                      "Lower Crooked", "Upper Crooked", "Little Deschutes", "Beaver-South Fork",
                                                                      "North Fork John Day", "Middle Fork John Day", "Upper John Day", "Klickitat") &
                                        !(event_site_name %in% origin_sites) &
                                        # Make sure to remove all instream and dams
                                        !grepl("Fishways", event_site_name) &
                                        !grepl("McNary Dam", event_site_name) &
                                        !grepl("Dalles Dam", event_site_name) &
                                        !grepl("Columbia River", event_site_name) &
                                        !grepl("John Day Dam", event_site_name) &
                                        !grepl("JDA", event_site_name) &
                                        !grepl("JO1|JO2", event_site_name) |
                                        event_site_name == "ROCK2C - Rock Creek, Columbia River (WA)")

BON_MCN_other_trib_sites <- BON_MCN_other_trib_sites_df$event_site_name


# Straying sites between MCN and branches (PRA or ICH):
# I actually don't think there are any - since the only subbasins are the Yakima and Walla Walla, and we care about those
# Basins: Yakima
# Subbasins: Walla Walla

# Straying sites between ICH and LGR:
# Basins: NA
# Subbasins: Palouse (but there aren't even any)
# Tucannon is already in the origins we are interested in
ICH_LGR_other_trib_sites_df <- subset(complete_event_site_metadata, 
                                      event_site_subbasin_name %in% c("Palouse") |
                                        event_site_name %in% c("ALMOTC - Almota Creek - tributary to Snake River",
                                                               "PENAWC - Penawawa Creek - tributary to Snake River"))

ICH_LGR_other_trib_sites <- ICH_LGR_other_trib_sites_df$event_site_name

# Straying sites upstream of LGR:
# Basins: Salmon
# Subbasins: Upper Grande Ronde, Clearwater, Lower Grande Ronde
# It looks like every one of these is contained within our tributaries of interest...

LGR_upstream_other_trib_sites_df <- subset(complete_event_site_metadata, 
                                           event_site_subbasin_name %in% c("Upper Grande Ronde", "Lower Grande Ronde", "Clearwater", "Lower North Fork Clearwater",
                                                                           "Upper North Fork Clearwater", "Lochsa", "Lower Selway", "Upper Selway", "Lower Salmon",
                                                                           "Little Salmon", "Middle Salmon-Chamberlain", "Lower Middle Fork Salmon", 
                                                                           "Upper Middle Fork Salmon", "Middle Salmon-Panther", "Pahsimeroi", "Lemhi") & 
                                             !(event_site_name %in% origin_sites) |
                                             event_site_name %in% c("ALPOWC - Alpowa Creek, lower Snake River",
                                                                    "TENMC2 - Tenmile Creek, tributary to Snake River"
                                             ) &
                                             !(event_site_name %in% origin_sites))

LGR_upstream_other_trib_sites <- LGR_upstream_other_trib_sites_df$event_site_name

# Straying sites between PRA and RIS
# Looking at the map, I don't think that there are any


# Straying sites between RIS and RRE
# I also don't think there are any here - there's only the Wenatchee, and we're interested in that


# Straying sites between RRE and WEL
# Again, nothing except the Entiat, which we care about

# Straying sites upstream of WEL
WEL_upstream_other_trib_sites_df <- subset(complete_event_site_metadata, 
                                           grepl("Nespelem", event_site_name) &
                                             !(event_site_name %in% origin_sites) |
                                             event_site_subbasin_name %in% c("Sanpoil", "Lower Spokane", "Franklin D. Roosevelt Lake")
                                           &
                                             !(event_site_name %in% origin_sites) |
                                             event_site_name %in% c("FST - Foster Creek"))

WEL_upstream_other_trib_sites <- WEL_upstream_other_trib_sites_df$event_site_name



##### In-river detection sites #####

pre_BON_inriver <- c("ESANIS - East Sand Island, Columbia River", "TWX - Estuary Towed Array (Exp.)",
                     "BONH - Bonneville Hatchery","PD7 - Columbia River Estuary rkm 70")
# This includes some dams that we are ignoring
BON_MCN_inriver <- c("COLR4 - Columbia River - Bonneville Dam to John Day Dam (km 234-347)",
                     "The Dalles Adult Fishways (combined)", "JDJ - John Day Dam Juvenile",
                     "LMILIS - Little Miller Island, Columbia River", "TDLPI - Lone Pine Island and associated unnamed islands near The Dalles Dam",
                     "COLR5 - Columbia River - John Day Dam to Snake River (km 347-522)",
                     # These are all new since 2017
                     "JO2 - John Day North Fish Ladder", "JDALD1 - JDA - Release into south fish ladder", 
                     "JO1 - John Day South Fish Ladder", "JDALD2 - JDA - Release into north fish ladder",
                     "John Day Dam Adult Fishways (combined)")
MCN_ICH_PRA_inriver <- c("BADGEI - Badger Island, Columbia River", "PRH - Priest Rapids Hatchery Outfall",
                         "RSH - Ringold Springs Hatch. Outfall")
# This includes dams (LMO and LGO) that we are ignoring
# These are some fallback routes - but because we aren't looking at these
# dams specifically, we can consider them in river detection sites
ICH_LGR_inriver <- c("GRJ - Lower Granite Dam Juvenile", "GOJ - Little Goose Dam Juvenile",
                     "LMJ - Lower Monumental Dam Juvenile",
                     "LGRTAL - LGR - Release into the Tailrace within 0.5 km downstream of Dam",
                     "LMA - Lower Monumental Adult Ladders", "GOA - Little Goose Fish Ladder",
                     "LYFE - Lyons Ferry Hatchery", "SNAKE1 - Snake River - mouth to Palouse River (km 0-96)",
                     "SNAKE2 - Snake River - Palouse River to Clearwater River (km 96-224)")

RRE_WEL_inriver <- c("WELH - Wells Hatchery", "WEH - Wells Dam Hatchery", "WELTAL - WEL - Release into the Tailrace within 0.5 km downstream of Dam",
                     "EBO - East Bank Hatchery Outfall")

upstream_WEL_inriver <- c("CHJO - Chief Joseph Hatchery", "OXBO - Oxbow Hatchery (IDFG)",
                          "COLR8 - Columbia River - Chelan Falls, WA to Grand Coulee Dam (km 809-960)",
                          "WELFBY - WEL - Release into the Forebay within 0.5 km upstream of Dam")

upstream_LGR_inriver <- c("SNAKE4 - Snake River - Salmon River to Hells Canyon Dam (km 303-397)")

##### Sites that are confirmed fallback (only certain dams) ####
BON_fallback_arrays <- c("BCC - BON PH2 Corner Collector", "B2J - Bonneville PH2 Juvenile")
MCN_fallback_arrays <- c("MCJ - McNary Dam Juvenile")
LGR_fallback_arrays <- c("GRJ - Lower Granite Dam Juvenile", "GRS - Lower Granite Dam Spillway")
RRE_fallback_arrays <- c("RRJ - Rocky Reach Dam Juvenile")
RIS_fallback_arrays <- c("RI2BYP - RIS - Release into the PH2 Juvenile Facility Bypass Flume/Pipe")
WEL_fallback_arrays <- c("WEJ - Wells Dam Bypass Bay Sample")

# Confirm that all sites have been categorized
setdiff(complete_event_site_metadata$event_site_name, 
        c(columbia_dams, snake_dams, # Adult fishways at dams
          origin_sites, # Natal origins
          LGR_upstream_other_trib_sites, ICH_LGR_other_trib_sites, # Other tributary sitses
          BON_MCN_other_trib_sites, WEL_upstream_other_trib_sites, # Other tributary sites continued
          pre_BON_inriver, BON_MCN_inriver, MCN_ICH_PRA_inriver, ICH_LGR_inriver, 
          upstream_WEL_inriver, RRE_WEL_inriver, upstream_LGR_inriver, #in river arrays
          BON_fallback_arrays, MCN_fallback_arrays, LGR_fallback_arrays, 
          RIS_fallback_arrays, RRE_fallback_arrays, WEL_fallback_arrays)) # fallback arrays

# Turn this into a dataframe (this was a silly way to do this)
site_classification <- data.frame(event_site_name = c(columbia_dams, snake_dams, # Adult fishways at dams
                                                      ASO_sites,  # natal origins
                                                      CLE_sites,
                                                      DES_sites,
                                                      ENT_sites,
                                                      FIF_sites,
                                                      GRRO_sites,
                                                      HOOD_sites,
                                                      IMN_sites,
                                                      JDR_sites,
                                                      MET_sites,
                                                      OKA_sites,
                                                      SAL_sites,
                                                      TUC_sites,
                                                      UMA_sites,
                                                      WAWA_sites,
                                                      WEN_sites,
                                                      YAK_sites,
                                                      LGR_upstream_other_trib_sites, ICH_LGR_other_trib_sites, # other trib sites
                                                      BON_MCN_other_trib_sites, WEL_upstream_other_trib_sites,
                                                      pre_BON_inriver, BON_MCN_inriver, MCN_ICH_PRA_inriver, ICH_LGR_inriver, #in river arrays
                                                      upstream_WEL_inriver, upstream_LGR_inriver, RRE_WEL_inriver,
                                                      BON_fallback_arrays, MCN_fallback_arrays, LGR_fallback_arrays, # fallback arrays
                                                      RRE_fallback_arrays, RIS_fallback_arrays, WEL_fallback_arrays,
                                                      "lost"),
                                  site_class = c(c("BON (adult)", "MCN (adult)", "PRA (adult)", "RIS (adult)",
                                                   "RRE (adult)", "WEL (adult)"),
                                                 c("ICH (adult)", "LGR (adult)"), 
                                                 rep("ASO_sites", length(ASO_sites)),
                                                 rep("CLE_sites", length(CLE_sites)),
                                                 rep("DES_sites", length(DES_sites)),
                                                 rep("ENT_sites", length(ENT_sites)),
                                                 rep("FIF_sites", length(FIF_sites)),
                                                 rep("GRRO_sites", length(GRRO_sites)),
                                                 rep("HOOD_sites", length(HOOD_sites)),
                                                 rep("IMN_sites", length(IMN_sites)),
                                                 rep("JDR_sites", length(JDR_sites)),
                                                 rep("MET_sites", length(MET_sites)),
                                                 rep("OKA_sites", length(OKA_sites)),
                                                 rep("SAL_sites", length(SAL_sites)),
                                                 rep("TUC_sites", length(TUC_sites)),
                                                 rep("UMA_sites", length(UMA_sites)),
                                                 rep("WAWA_sites", length(WAWA_sites)),
                                                 rep("WEN_sites", length(WEN_sites)),
                                                 rep("YAK_sites", length(YAK_sites)),
                                                 rep("LGR_upstream_other_trib_sites", length(LGR_upstream_other_trib_sites)),
                                                 rep("ICH_LGR_other_trib_sites", length(ICH_LGR_other_trib_sites)),
                                                 rep("BON_MCN_other_trib_sites", length(BON_MCN_other_trib_sites)),
                                                 rep("WEL_upstream_other_trib_sites", length(WEL_upstream_other_trib_sites)),
                                                 rep("pre_BON_inriver", length(pre_BON_inriver)), 
                                                 rep("BON_MCN_inriver", length(BON_MCN_inriver)),
                                                 rep("MCN_ICH_PRA_inriver", length(MCN_ICH_PRA_inriver)),
                                                 rep("ICH_LGR_inriver", length(ICH_LGR_inriver)), 
                                                 rep("upstream_WEL_inriver", length(upstream_WEL_inriver)), 
                                                 rep("upstream_LGR_inriver", length(upstream_LGR_inriver)), 
                                                 rep("RRE_WEL_inriver", length(RRE_WEL_inriver)), 
                                                 rep("BON_fallback_arrays", length(BON_fallback_arrays)),
                                                 rep("MCN_fallback_arrays", length(MCN_fallback_arrays)),
                                                 rep("LGR_fallback_arrays", length(LGR_fallback_arrays)),
                                                 rep("RRE_fallback_arrays", length(RRE_fallback_arrays)),
                                                 rep("RIS_fallback_arrays", length(RIS_fallback_arrays)),
                                                 rep("WEL_fallback_arrays", length(WEL_fallback_arrays)),
                                                 "lost"))

# Create a new dataframe for state (location in the system)
site_classification %>%
  mutate(
    state = ifelse(
      site_class == "BON (adult)", "mainstem, BON to MCN",
      ifelse(
        site_class == "MCN (adult)", "mainstem, MCN to ICH or PRA",
        ifelse(
          site_class == "PRA (adult)", "mainstem, PRA to RIS",
          ifelse(
            site_class == "RIS (adult)", "mainstem, RIS to RRE",
            ifelse(
              site_class == "RRE (adult)", "mainstem, RRE to WEL",
              ifelse(
                site_class == "WEL (adult)", "mainstem, upstream of WEL",
                ifelse(
                  site_class == "ICH (adult)", "mainstem, ICH to LGR",
                  ifelse(
                    site_class == "LGR (adult)", "mainstem, upstream of LGR",
                    ifelse(
                      # Natal origins
                      site_class == "ASO_sites", "Asotin Creek",
                      ifelse(
                        site_class == "CLE_sites", "Clearwater River",
                        ifelse(
                          site_class == "DES_sites", "Deschutes River",
                          ifelse(
                            site_class == "ENT_sites", "Entiat River",
                            ifelse(
                              site_class == "FIF_sites", "Fifteenmile Creek",
                              ifelse(
                                site_class == "GRRO_sites", "Grande Ronde River",
                                ifelse(
                                  site_class == "HOOD_sites", "Hood River",
                                  ifelse(
                                  site_class == "IMN_sites", "Imnaha River",
                                  ifelse(
                                    site_class == "JDR_sites", "John Day River",
                                    ifelse(
                                      site_class == "MET_sites", "Methow River",
                                      ifelse(
                                        site_class == "OKA_sites", "Okanogan River",
                                        ifelse(
                                          site_class == "SAL_sites", "Salmon River",
                                          ifelse(
                                            site_class == "TUC_sites", "Tucannon River",
                                            ifelse(
                                              site_class == "UMA_sites", "Umatilla River",
                                              ifelse(
                                                site_class == "WAWA_sites", "Walla Walla River",
                                                ifelse(
                                                  site_class == "WEN_sites", "Wenatchee River",
                                                  ifelse(
                                                    site_class == "YAK_sites", "Yakima River",
                                                    ifelse(
                                                      
                                                      # Other tributaries
                                                      site_class == "LGR_upstream_other_trib_sites", "Upstream LGR other tributaries",
                                                      ifelse(
                                                        site_class == "ICH_LGR_other_trib_sites", "ICH to LGR other tributaries",
                                                        ifelse(
                                                          site_class == "BON_MCN_other_trib_sites", "BON to MCN other tributaries",
                                                          ifelse(
                                                            site_class == "WEL_upstream_other_trib_sites", "Upstream WEL other tributaries",
                                                            ifelse(
                                                              
                                                              # In river sites
                                                              site_class == "pre_BON_inriver", "mainstem, mouth to BON",
                                                              ifelse(
                                                                site_class == "BON_MCN_inriver", "mainstem, BON to MCN",
                                                                ifelse(
                                                                  site_class == "MCN_ICH_PRA_inriver", "mainstem, MCN to ICH or PRA",
                                                                  ifelse(
                                                                    site_class == "ICH_LGR_inriver", "mainstem, ICH to LGR",
                                                                    ifelse(
                                                                      site_class == "upstream_WEL_inriver", "mainstem, upstream of WEL",
                                                                      ifelse(
                                                                        site_class == "upstream_LGR_inriver", "mainstem, upstream of LGR",
                                                                        ifelse(
                                                                          site_class == "RRE_WEL_inriver", "mainstem, RRE to WEL", 
                                                                          ifelse(
                                                                          site_class == "BON_fallback_arrays", "mainstem, mouth to BON",
                                                                          ifelse(
                                                                            site_class == "MCN_fallback_arrays", "mainstem, BON to MCN",
                                                                            ifelse(
                                                                              site_class == "LGR_fallback_arrays", "mainstem, ICH to LGR", 
                                                                              ifelse(
                                                                                site_class == "RRE_fallback_arrays", "mainstem, RIS to RRE", 
                                                                                ifelse(
                                                                                  site_class == "RIS_fallback_arrays", "mainstem, PRA to RIS", 
                                                                                  ifelse(
                                                                                    site_class == "WEL_fallback_arrays", "mainstem, RRE to WEL", 
                                                                                    ifelse(site_class == "lost", "lost", NA)
                                                                                  ))))))))))))))))))))))))))))))))))))))))))) -> site_classification

# Export site classification
write.csv(site_classification, "site_classification.csv")
# write.csv(site_classification, here::here("covariate_data", "site_classification.csv"))


# Add info to detection history
complete_det_hist <- read.csv("complete_det_hist.csv", row.names = 1)

complete_det_hist %>% 
  left_join(., site_classification, by = "event_site_name") -> det_hist


# Turn this into state history, not detection history.
# This needs to include that is implicit in the detection history; i.e., if
# a fish was seen in the BON adult ladders in consecutive detection events,
# it must have fallen back over and thus entered the pre BON mainstem state

# So this code currently works to detect some implicit fallback, i.e., detection
# at the same dam twice in a row. However, it isn't yet able to capture implicit
# fallback when it is detected in a tributary. For example, if a fish is seen
# at McNary, and then in the John Day, it doesn't register that it was in the
# mainstem between BON and MCN between those times. But perhaps we can
# change that after the fact?
# Between two dams, it should always be BON -> mainstem -> (tributary) -> mainstem -> MCN
# Don't have to go to tributary, but if they do they must go to mainstem before or after

# Need a function that can insert rows into the DF for the implicit states

# Order of sites (no tributaries), Columbia River
site_order_notrib_columbia <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                "mainstem, upstream of WEL")

# Order of sites (no tributaries), Snake
# site_order_notrib_snake <- c("mainstem, ICH to LGR", "mainstem, upstream of LGR")
# Contains all non-tributary sites from upstream of LGR to downstream of BON
site_order_notrib_snake <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                             "mainstem, MCN to ICH or PRA",
                             "mainstem, ICH to LGR",
                             "mainstem, upstream of LGR")


# Site order, no tributaries, Snake directly to upper Columbia
site_order_notrib_columbia_snake <- c("mainstem, upstream of WEL",
                                      "mainstem, RRE to WEL", 
                                      "mainstem, RIS to RRE",
                                      "mainstem, PRA to RIS",
                                      "mainstem, MCN to ICH or PRA",
                                      "mainstem, ICH to LGR",
                                      "mainstem, upstream of LGR")


# Upper Columbia sites
upper_columbia_sites <- c("mainstem, PRA to RIS",
                          "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                          "mainstem, upstream of WEL")

upper_columbia_tribs <- c("Okanogan River", "Methow River", "Upstream WEL other tributaries", "Entiat River", "Wenatchee River")
# Make sure they're right
# setdiff(upper_columbia_tribs, tributary_mainstem$tributary)

# Snake River sites
# snake_sites <- c("mainstem, ICH to LGR", 
#                             "ICH to LGR tributaries", "mainstem, ICH to LGR",
#                             "mainstem, upstream of LGR", "Upstream LGR tributaries",
#                             "mainstem, upstream of LGR")
snake_sites <- c("mainstem, ICH to LGR", 
                 "mainstem, upstream of LGR")
snake_tribs <- c("Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River",
                 "Tucannon River", "Upstream LGR other tributaries", "ICH to LGR other tributaries")
# setdiff(snake_tribs, tributary_mainstem$tributary)


# Shared sites (lower Columbia, before confluence of Snake and Columbia)
# shared_sites <- c("mainstem, mouth to BON", "mainstem, BON to MCN",
#                                      "natal tributaries",
#                                      "BON to MCN tributaries",
#                   "MCN to PRA or ICH tributaries",
#                                      "mainstem, MCN to ICH or PRA")
shared_sites <- c("mainstem, mouth to BON", "mainstem, BON to MCN",
                  "mainstem, MCN to ICH or PRA")

shared_tribs <- c("Deschutes River", "BON to MCN other tributaries", "John Day River",
                  "Hood River", "Fifteenmile Creek", "Umatilla River", "Yakima River",
                  "Walla Walla River")
# setdiff(shared_tribs, tributary_mainstem$tributary)

# Make sure we're capturing all tribs
# setdiff(tributary_mainstem$tributary, c(upper_columbia_tribs, snake_tribs, shared_tribs))

# Tributary df - tributaries and what part of the mainstem they're in
# We will use this to determine where they need to move to next, and then
# use the no tributary site order to fill in the rest of the history

tributary_mainstem <- data.frame(tributary = c("Asotin Creek", 
                                               "Clearwater River",
                                               "Deschutes River", 
                                               "Entiat River", 
                                               "Fifteenmile Creek", 
                                               "Grande Ronde River", 
                                               "Hood River",
                                               "Imnaha River",
                                               "John Day River", 
                                               "Methow River", 
                                               "Okanogan River", 
                                               "Salmon River", 
                                               "Tucannon River", 
                                               "Umatilla River",
                                               "Walla Walla River",
                                               "Wenatchee River", 
                                               "Yakima River",
                                               "Upstream LGR other tributaries",
                                               "ICH to LGR other tributaries",
                                               "BON to MCN other tributaries",
                                               "Upstream WEL other tributaries"
),
mainstem = c("mainstem, upstream of LGR",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, RRE to WEL",
             "mainstem, BON to MCN",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of WEL",
             "mainstem, upstream of WEL",
             "mainstem, upstream of LGR",
             "mainstem, ICH to LGR",
             "mainstem, BON to MCN",
             "mainstem, MCN to ICH or PRA",
             "mainstem, RIS to RRE",
             "mainstem, MCN to ICH or PRA",
             "mainstem, upstream of LGR",
             "mainstem, ICH to LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of WEL"))


# Inspect detection histories after tributary detections
# Do we want to remove this? Probably indicates post-spawning behavior
# det_hist %>% 
#   group_by(tag_code) %>% 
#   subset(state %in% tributary_mainstem$tributary) %>% 
#   filter(end_time == max(end_time)) %>% 
#   dplyr::select(tag_code, end_time) %>% 
#   dplyr::rename(spawning = end_time) -> spawning_times

# complete_det_hist %>% 
#   left_join(., spawning_times, by = "tag_code") %>%  
#   group_by(tag_code) %>% 
#   filter(end_time > spawning) -> post_spawn_history

# A lot of these do seem like iteroparous individuals migrating back downstream. But not all of them.
# Would bias our fallback estimates to include them! But we need to be flexible enough to 
# include this behavior, since some of it is tributary dip-ins

##### FOR LOOP #####

# Make a function to insert rows for implicit detections
# This works for the det_hist dataframe
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}


# Create empty df to store states
stepwise_states <- data.frame(tag_code = character(),
                              state = character(),
                              date_time = as.POSIXct(character()),
                              pathway = character())


# Before running this, save the original detection history
# det_hist_original <- det_hist
# Reset detection history
# det_hist <- det_hist_original


# Add a counter to keep track of the number of added rows.
# This isn't actually needed for the code after all, but is still interesting
added_rows <- 0

# Our for loop is too long, and will give us an error. However, we need to
# make sure that we are accounting for the rows that we are inserting into
# our dataframe during this for loop, so to be safe we are doubling the 
# length of the for loop. It will run fine, it will just give an error once
# the dataframe ends, but that doesn't matter for our purposes.
# NOTE: If you don't get an error, that means your for loop isn't long enough
# and you need to extend it
# The error should read: 
# "Error in if (det_hist[i, "tag_code"] == det_hist[i - 1, "tag_code"]) { : 
# missing value where TRUE/FALSE needed

# Start a counter for which tag code
fish_counter <- 1

nfish <- length(unique(det_hist$tag_code))

# Append a dummy fish

det_hist %>% 
  bind_rows(data.frame(tag_code = "dummy_fish",
              event_site_name = NA,
              start_time = NA,
              end_time = NA,
              event_site_basin_name = NA,
              event_site_subbasin_name = NA,
              event_site_latitude = NA,
              event_site_longitude = NA,
              site_class = "dummy_site", 
              state = NA)) -> det_hist

# Add the "pathway" field
# det_hist %>% 
#   mutate(pathway = "NA") -> det_hist
  

# Loop through each row of the detection history
# Okay so for hyak I need to stop it before it gets an error because it reached the end of the df
# Use a dummy fish to stop it

# unique(det_hist$tag_code)[46848]

# Examine the detection histories for the funky fish

# for troubleshooting: load problem tags, subset

# Save the original
# det_hist_orig <- det_hist

# read.table("problem_tags.txt") %>% 
#   dplyr::rename("tag_code" = V1) -> problem_tags

# det_hist <- subset(det_hist_orig, tag_code %in% problem_tags$tag_code)

# reset the indices
rownames(det_hist) <- NULL

length(unique(det_hist$tag_code)) # 215 unique ones of 412 total

# Find individual fish:
# subset(det_hist, tag_code == "3DD.0077E5598A") %>% dplyr::select(tag_code, site_class, end_time)
# subset(states_complete, tag_code == "3DD.0077E5598A")


stepwise_states <- data.frame(tag_code = character(),
                              state = character(),
                              date_time = as.POSIXct(character()),
                              pathway = character())
# This tag code 384.1B796A520E is rows 780 - 795 of det hist
for (i in 1:round(1.5 * nrow(det_hist), 0)) {
# for (i in 88760:88800){
# for (i in 7200:7300) {
# for (i in 22500:23000){
# for (i in 31900:32500) {
# for (i in (nrow(det_hist)-20):(nrow(det_hist) + 20)) { # Test to see if dummy fish solution works
  print(paste0("Row #", i))
  if (fish_counter > nfish){
    # So here, if our fish counter is greater than the number of REAL fish in the dataset (this is due to our dummy fish) - do nothing
    print(paste0("Run complete - no more fish"))
  }
  # If it's the first entry, store it
  else if (i == 1) {
    # Store the tag code
    stepwise_states[i, 'tag_code'] <-
      det_hist[i, 'tag_code']
    # Store the current state
    stepwise_states[i, 'state'] <- det_hist[i, 'state']
    # Store the time entering this state
    stepwise_states[i, 'date_time'] <-
      det_hist[i, 'end_time']
    # Store the transition site
    stepwise_states[i, 'pathway'] <-
      det_hist[i, 'site_class']
  }
  
  # If we've appended lots of entries and the stepwise states and det_hist dfs no longer line up, then skip
  # Try just that if it's an implicit entry, then skip it because it should already be in the stepwise states as well
  else if (det_hist[i, "site_class"] == "implicit"){
    # skip!
  }
  
  # For all other entries, look at the previous entry.
  # If it's the same fish as the previous entry:
  else if (det_hist[i, 'tag_code'] == det_hist[i - 1, 'tag_code']) {
    # If it's in a different state than it was previously, record it
    if (det_hist[i, 'state'] != det_hist[i - 1, 'state']) {
      ### NEW CONDITIONS: COMPLETE DETECTION HISTORY (NOT REPEAT DAM DET)
      # We have to account for what happens if we see a fish at an
      # site, and then next at a different site while
      # skipping the sites in between (could be upstream or downstream).
      
      
      # if it's in the right part of the mainstem before, then we skip it.
      # else (anything else) means it's not. In that case, we have to determine
      # what part of the mainstem it needs to have been in to get to that tributary,
      # then use that information to fill in the missing sites from the order of
      # sites DFs
      
      ### TRIBUTARIES
      # If it's in a tributary and wasn't in the right part of the mainstem prev:
      # Also, must not already be in the tributary
      if (det_hist[i, 'state'] %in% tributary_mainstem$tributary) {
        
        # SPECIAL CASE: IF THERE WERE CONSECUTIVE TRIBUTARY DETECTIONS
        
        # Edits 07-16-2022: I think that the code within the if statements is correct, the if statement itself just isn't right.
        # Changes: changed the if statements to look for tributaries by part of river
        
        if (det_hist[i, 'state'] %in% tributary_mainstem$tributary &
            det_hist[i - 1, 'state'] %in% tributary_mainstem$tributary){
          
          # Zeroth: Jumping between Snake and Columbia, without visiting
          # shared sites in between. We see this with some fish, if they are
          # seen at PRA and then ICH or vice versa, without anything in between
          if (det_hist[i - 1, 'state'] %in% upper_columbia_tribs &
              det_hist[i, 'state'] %in% snake_tribs |
              det_hist[i - 1, 'state'] %in% snake_tribs &
              det_hist[i, 'state'] %in% upper_columbia_tribs) {
            
            ####### THIS SECTION IDENTICAL, JUST WITH SHARED, COLUMBIA, and SNAKE
            # Get the current mainstem site & index
            mainstem_site_current <- subset(tributary_mainstem, tributary == det_hist[i, 'state'])$mainstem
            mainstem_index_current <- which(site_order_notrib_columbia_snake %in% mainstem_site_current)
            
            # Get the previous mainstem site & index
            mainstem_site_previous <- subset(tributary_mainstem, tributary == det_hist[i-1, 'state'])$mainstem
            mainstem_index_previous <- which(site_order_notrib_columbia_snake %in% mainstem_site_previous)
            
            # If the fish is moving upstream
            if (mainstem_index_current > mainstem_index_previous){
              
              # sequence from current to previous index
              index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
              
              # Get the missing sites
              missing_sites <- site_order_notrib_columbia_snake[index_order]
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order))) {
                ### EDIT 2022-07-16
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
                missing_sites <- site_order_notrib_columbia_snake[index_order]
                ### end edit
                
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <- missing_sites[j]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                index_order <- seq(mainstem_index_current, mainstem_index_previous, by = -1)
                
                missing_sites <- site_order_notrib_columbia_snake[index_order]
                
                missing_site <- missing_sites[j]
                ### End edit
                
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
            # If the fish is moving downstream
            if (mainstem_index_current < mainstem_index_previous){
              # sequence from current to previous index
              index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
              
              # Get the missing sites
              missing_sites <- site_order_notrib_columbia_snake[index_order]
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order))) {
                ### EDIT 2022-07-16
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
                missing_sites <- site_order_notrib_columbia_snake[index_order]
                ### end edit
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <- missing_sites[j]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                index_order <- seq(mainstem_index_current, mainstem_index_previous, by = 1)
                
                missing_sites <- site_order_notrib_columbia_snake[index_order]
                
                missing_site <- missing_sites[j]
                ### End edit
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
          }
          
          
          
          
          
          
          
          
          # First: Shared sites - need two options, depending on where it previously was
          # Shared sites can use either Columbia or Snake sites
          # Condition: If the indices are off by more than one (plus or minus one)
          else if (det_hist[i, 'state'] %in% shared_tribs) {
            
            # EDIT 2022-07-16: adding the if statement for if in columbia or snake
            # If it was previously in the Columbia
            if (det_hist[i - 1, 'state'] %in% c(shared_tribs, upper_columbia_tribs)){
              ####### THIS SECTION IDENTICAL, JUST WITH SHARED, COLUMBIA, and SNAKE
              # Get the current mainstem site & index
              mainstem_site_current <- subset(tributary_mainstem, tributary == det_hist[i, 'state'])$mainstem
              mainstem_index_current <- which(site_order_notrib_columbia %in% mainstem_site_current)
              
              # Get the previous mainstem site & index
              mainstem_site_previous <- subset(tributary_mainstem, tributary == det_hist[i-1, 'state'])$mainstem
              mainstem_index_previous <- which(site_order_notrib_columbia %in% mainstem_site_previous)
              
              # If the fish is moving upstream
              if (mainstem_index_current > mainstem_index_previous){
                
                # sequence from current to previous index
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
                
                # Get the missing sites
                missing_sites <- site_order_notrib_columbia[index_order]
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order))) {
                  ### EDIT 2022-07-16
                  index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
                  missing_sites <- site_order_notrib_columbia[index_order]
                  ### end edit
                  
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <- missing_sites[j]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                  index_order <- seq(mainstem_index_current, mainstem_index_previous, by = -1)
                  
                  missing_sites <- site_order_notrib_columbia[index_order]
                  
                  missing_site <- missing_sites[j]
                  ### End edit
                  
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
                
              }
              
              # If the fish is moving downstream
              if (mainstem_index_current < mainstem_index_previous){
                # sequence from current to previous index
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
                
                # Get the missing sites
                missing_sites <- site_order_notrib_columbia[index_order]
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order))) {
                  ### EDIT 2022-07-16
                  index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
                  missing_sites <- site_order_notrib_columbia[index_order]
                  ### end edit
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <- missing_sites[j]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                  index_order <- seq(mainstem_index_current, mainstem_index_previous, by = 1)
                  
                  missing_sites <- site_order_notrib_columbia[index_order]
                  
                  missing_site <- missing_sites[j]
                  ### End edit
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
                
              }
              
            }
            
            # If it was previously in the Snake
            else{
              ####### THIS SECTION IDENTICAL, JUST WITH SHARED, COLUMBIA, and SNAKE
              # Get the current mainstem site & index
              mainstem_site_current <- subset(tributary_mainstem, tributary == det_hist[i, 'state'])$mainstem
              mainstem_index_current <- which(site_order_notrib_snake %in% mainstem_site_current)
              
              # Get the previous mainstem site & index
              mainstem_site_previous <- subset(tributary_mainstem, tributary == det_hist[i-1, 'state'])$mainstem
              mainstem_index_previous <- which(site_order_notrib_snake %in% mainstem_site_previous)
              
              # If the fish is moving upstream
              if (mainstem_index_current > mainstem_index_previous){
                
                # sequence from current to previous index
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
                
                # Get the missing sites
                missing_sites <- site_order_notrib_snake[index_order]
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order))) {
                  ### EDIT 2022-07-16
                  index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
                  missing_sites <- site_order_notrib_snake[index_order]
                  ### end edit
                  
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <- missing_sites[j]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                  index_order <- seq(mainstem_index_current, mainstem_index_previous, by = -1)
                  
                  missing_sites <- site_order_notrib_snake[index_order]
                  
                  missing_site <- missing_sites[j]
                  ### End edit
                  
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
                
              }
              
              # If the fish is moving downstream
              if (mainstem_index_current < mainstem_index_previous){
                # sequence from current to previous index
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
                
                # Get the missing sites
                missing_sites <- site_order_notrib_snake[index_order]
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order))) {
                  ### EDIT 2022-07-16
                  index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
                  missing_sites <- site_order_notrib_snake[index_order]
                  ### end edit
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <- missing_sites[j]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                  index_order <- seq(mainstem_index_current, mainstem_index_previous, by = 1)
                  
                  missing_sites <- site_order_notrib_snake[index_order]
                  
                  missing_site <- missing_sites[j]
                  ### End edit
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
                
              }
            }

            
          }
          
          
          
          
          
          
          
          # Second: Upper Columbia sites
          else if (det_hist[i, 'state'] %in% upper_columbia_tribs) {
            ####### THIS SECTION IDENTICAL, JUST WITH SHARED, COLUMBIA, and SNAKE
            # Get the current mainstem site & index
            mainstem_site_current <- subset(tributary_mainstem, tributary == det_hist[i, 'state'])$mainstem
            mainstem_index_current <- which(site_order_notrib_columbia %in% mainstem_site_current)
            
            # Get the previous mainstem site & index
            mainstem_site_previous <- subset(tributary_mainstem, tributary == det_hist[i-1, 'state'])$mainstem
            mainstem_index_previous <- which(site_order_notrib_columbia %in% mainstem_site_previous)
            
            # If the fish is moving upstream
            if (mainstem_index_current > mainstem_index_previous){
              
              # sequence from current to previous index
              index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
              
              # Get the missing sites
              missing_sites <- site_order_notrib_columbia[index_order]
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order))) {
                ### EDIT 2022-07-16
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
                missing_sites <- site_order_notrib_columbia[index_order]
                ### end edit
                
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <- missing_sites[j]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                index_order <- seq(mainstem_index_current, mainstem_index_previous, by = -1)
                
                missing_sites <- site_order_notrib_columbia[index_order]
                
                missing_site <- missing_sites[j]
                ### End edit
                
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                

                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
            # If the fish is moving downstream
            if (mainstem_index_current < mainstem_index_previous){
              # sequence from current to previous index
              index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
              
              # Get the missing sites
              missing_sites <- site_order_notrib_columbia[index_order]
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order))) {
                ### EDIT 2022-07-16
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
                missing_sites <- site_order_notrib_columbia[index_order]
                ### end edit
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <- missing_sites[j]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                index_order <- seq(mainstem_index_current, mainstem_index_previous, by = 1)
                
                missing_sites <- site_order_notrib_columbia[index_order]
                
                missing_site <- missing_sites[j]
                ### End edit
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
          }
          
          
          
          # Third: Snake River sites
          else if (det_hist[i, 'state'] %in% snake_tribs) {
            ####### THIS SECTION IDENTICAL, JUST WITH SHARED, COLUMBIA, and SNAKE
            # Get the current mainstem site & index
            mainstem_site_current <- subset(tributary_mainstem, tributary == det_hist[i, 'state'])$mainstem
            mainstem_index_current <- which(site_order_notrib_snake %in% mainstem_site_current)
            
            # Get the previous mainstem site & index
            mainstem_site_previous <- subset(tributary_mainstem, tributary == det_hist[i-1, 'state'])$mainstem
            mainstem_index_previous <- which(site_order_notrib_snake %in% mainstem_site_previous)
            
            # If the fish is moving upstream
            if (mainstem_index_current > mainstem_index_previous){
              
              # sequence from current to previous index
              index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
              
              # Get the missing sites
              missing_sites <- site_order_notrib_snake[index_order]
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order))) {
                ### EDIT 2022-07-16
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = 1)
                missing_sites <- site_order_notrib_snake[index_order]
                ### end edit
                
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <- missing_sites[j]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                index_order <- seq(mainstem_index_current, mainstem_index_previous, by = -1)
                
                missing_sites <- site_order_notrib_snake[index_order]
                
                missing_site <- missing_sites[j]
                ### End edit
                
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
            # If the fish is moving downstream
            if (mainstem_index_current < mainstem_index_previous){
              # sequence from current to previous index
              index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
              
              # Get the missing sites
              missing_sites <- site_order_notrib_snake[index_order]
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order))) {
                ### EDIT 2022-07-16
                index_order <- seq(mainstem_index_previous, mainstem_index_current, by = -1)
                missing_sites <- site_order_notrib_snake[index_order]
                ### end edit
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <- missing_sites[j]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                ### EDIT 2022-07-16: flipping the order of site here (as we have done in other parts of the for loop)
                index_order <- seq(mainstem_index_current, mainstem_index_previous, by = 1)
                
                missing_sites <- site_order_notrib_snake[index_order]
                
                missing_site <- missing_sites[j]
                ### End edit
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
          }
          
          
          
          
        
      
        
      }
        
        
        # Now if it wasn't consecutive tributary detections, make sure the mainstem detections make sense
      
        else if (det_hist[i, 'state'] %in% tributary_mainstem$tributary &
            (det_hist[i - 1, 'state'] !=  subset(tributary_mainstem, tributary ==
                                                  det_hist[i, 'state'])$mainstem
            )) {
          mainstem_site <-
            subset(tributary_mainstem, tributary == det_hist[i, 'state'])$mainstem
          
          # INSERT ROW FOR THE MAINSTEM SITE
          # Insert a new row into stepwise states, with the implicit detection site
          # Tag code, state, and time (which is NA)
          implicit_state_1 <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                         state = mainstem_site,
                                         date_time = NA,
                                         pathway = "implicit")
          #1.26.22#
          # Need to append zeros, but only if i > nrow
          if (i > nrow(stepwise_states)){
            NA_rows <- data.frame(tag_code = rep(as.character(NA), i - (nrow(stepwise_states) + 1)),
                                  state = rep(as.character(NA), i - (nrow(stepwise_states) + 1)),
                                  date_time = rep(as.POSIXct(as.character(NA)), i - (nrow(stepwise_states) + 1)),
                                  pathway = rep(as.character(NA), i - (nrow(stepwise_states) + 1)))
            
            stepwise_states %>%
              # Add rows of NAs for any rows that were the same state
              bind_rows(., NA_rows) %>% 
              bind_rows(., implicit_state_1) -> stepwise_states
          }
          
          else{
            stepwise_states %>%
              bind_rows(., implicit_state_1) -> stepwise_states
          }
          # # # 
          
          # stepwise_states <-
          #   insertRow(stepwise_states, implicit_state_1, i)
          
          # Insert a new row into original detection history, with
          # implicit detection site info
          implicit_detection <- c(det_hist[i, 'tag_code'],
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  "implicit",
                                  mainstem_site)
          
          # Also, change the value in the original dataframe
          det_hist <-
            insertRow(det_hist, implicit_detection, i)
          
          # Add 1 to the number of added rows
          added_rows <- added_rows + 1
          
          
          # Add lines to detect missing sites, in shared, Columbia, or Snake
          ################################################################################################
          # Zeroth: Jumping between Snake and Columbia, without visiting
          # shared sites in between. We see this with some fish, if they are
          # seen at PRA and then ICH or vice versa, without anything in between
          if (det_hist[i - 1, 'state'] %in% upper_columbia_sites &
              det_hist[i, 'state'] %in% snake_sites |
              det_hist[i - 1, 'state'] %in% snake_sites &
              det_hist[i, 'state'] %in% upper_columbia_sites) {
            # These inherently have to skip, since they're not seen
            # at a shared site (MCN to ICH/PRA is the only route, and by
            # definition they're not seen there)
            # Get the index of the current site
            current_index <-
              which(site_order_notrib_columbia_snake %in% det_hist[i, 'state'])
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_columbia_snake %in% det_hist[i - 1, 'state'])
            
            
            # sequence from current to previous index
            if (current_index < previous_index) {
              # index_order <- seq(current_index, previous_index, by = 1)
              index_order <- seq(previous_index, current_index, by = -1)
              
            } else {
              index_order <- seq(previous_index, current_index, by = 1)
            }
            
            
            # Count the number of sites you need to add and loop through
            for (j in 1:(length(index_order) - 2)) {

              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_columbia_snake[index_order[1 + j]]
              implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                           state = missing_site,
                                           date_time = NA,
                                           pathway = "implicit")
              
              stepwise_states %>% 
                bind_rows(., implicit_state) -> stepwise_states
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              
              # Need to flip the order of sites for these - but it depends on the order of the sites
              if (current_index < previous_index){
                index_order <- seq(current_index, previous_index, by = 1)
              }
              else {
                index_order <- seq(current_index, previous_index, by = -1)
              }
              
              missing_site <- site_order_notrib_columbia_snake[index_order[j+1]]
              # missing_site <- site_order_notrib_columbia_snake[index_order[j]]
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              implicit_detection <- c(det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              det_hist <-
                insertRow(det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
            
            
          }
          
          
          # First: Shared sites - need two options, depending on where it previously was
          # Shared sites can use either Columbia or Snake sites
          # Condition: If the indices are off by more than one (plus or minus one)
          else if (det_hist[i, 'state'] %in% shared_sites) {
            # If it was previously in the Columbia
            if (det_hist[i - 1, 'state'] %in% c(shared_sites, upper_columbia_sites)) {
              # Make sure that if it was seen in the adult fishways at a dam, the
              # detection history contains the state downstream of that dam
              if (det_hist[i, 'site_class'] %in% c(
                "BON (adult)",
                "MCN (adult)",
                "PRA (adult)",
                "RIS (adult)",
                "RRE (adult)",
                "WEL (adult)",
                "ICH (adult)",
                "LGR (adult)"
              ) &
              (which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1) !=
              which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])) {
                
                # Get the missing index
                missing_index <-
                  which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
                
                ### EDIT 2022-07-16: moving this inside for loop
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                
                # Count the number of sites you need to add and loop through
                # for (j in 1:(length(index_order) - 2)) {
                # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                for (j in 1:(length(index_order) - 1)) {
                  ### EDIT 2022-7-16
                  # sequence from current to previous index
                  if (missing_index < previous_index) {
                    # index_order <- seq(current_index, previous_index, by = 1)
                    index_order <- seq(previous_index, missing_index, by = -1)
                    
                  } else {
                    index_order <- seq(previous_index, missing_index, by = 1)
                  }
                  # end edit
                  
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (missing_index < previous_index){
                    index_order <- seq(missing_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(missing_index, previous_index, by = -1)
                  }
                  
                  # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                  missing_site <- site_order_notrib_columbia[index_order[j]]
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
                
              }
              else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                       which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) + 1 |
                       which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                       which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) - 1) {
                # If we see it at two non-consecutive mainstem locations
                
                # If the next site skips sites, insert the missing sites
                # Get the index of the current site
                current_index <-
                  which(site_order_notrib_columbia %in% det_hist[i, 'state'])
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
                
                # sequence from current to previous index
                if (current_index < previous_index) {
                  index_order <- seq(current_index, previous_index, by = 1)
                  # index_order <- seq(previous_index, current_index, by = -1)
                  
                # } else {
                #   # index_order <- seq(previous_index, current_index, by = 1)
                #   index_order <- seq(previous_index, current_index, by = -1)
                # }
                } else {
                  # index_order <- seq(previous_index, current_index, by = 1)
                  index_order <- seq(current_index, previous_index, by = -1)
                }
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order) - 2)) {
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  ############
                  # stepwise_states %>% 
                  #   bind_rows(., implicit_state) -> stepwise_states
                  
                  # Also, change the value in the original dataframe
                  stepwise_states <-
                    insertRow(stepwise_states, implicit_state, i)
                  
                  ############
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  # if (missing_index < previous_index){
                  #   index_order <- seq(missing_index, previous_index, by = 1)
                  # }
                  # else {
                  #   index_order <- seq(missing_index, previous_index, by = -1)
                  # }
                  # 
                  # # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                  # missing_site <- site_order_notrib_columbia[index_order[j]]
                  
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
              }
              
              
              
            }
            
            # If it was previously in the Snake
            else if (det_hist[i - 1, 'state'] %in% snake_sites) {
              # Make sure that if it was seen in the adult fishways at a dam, the
              # detection history contains the state downstream of that dam
              if (det_hist[i, 'site_class'] %in% c(
                "BON (adult)",
                "MCN (adult)",
                "PRA (adult)",
                "RIS (adult)",
                "RRE (adult)",
                "WEL (adult)",
                "ICH (adult)",
                "LGR (adult)"
              ) &
              (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
              which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
                
                # Get the missing index
                missing_index <-
                  which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
                
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(missing_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                
                # Count the number of sites you need to add and loop through
                # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                for (j in 1:(length(index_order) - 1)) {
                  
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (missing_index < previous_index){
                    index_order <- seq(missing_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(missing_index, previous_index, by = -1)
                  }
                  
                  # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                  missing_site <- site_order_notrib_columbia[index_order[j]]
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
              }
              if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                  which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) + 1 |
                  which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                  which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) - 1) {
                # If we see it at two non-consecutive mainstem locations
                
                # If the next site skips sites, insert the missing sites
                # Get the index of the current site
                current_index <-
                  which(site_order_notrib_snake %in% det_hist[i, 'state'])
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
                
                ### EDIT 2022-07-16: Moving this inside the below for loop (same issue as before)
                # sequence from current to previous index
                if (current_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, current_index, by = -1)
                  
                # } else {
                #   index_order <- seq(previous_index, current_index, by = 1)
                # }
                } else {
                  # index_order <- seq(previous_index, current_index, by = 1)
                  index_order <- seq(current_index, previous_index, by = -1)
                }
                
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order) - 2)) {
                  ### EDIT 2022-07-16
                  if (current_index < previous_index) {
                    # index_order <- seq(current_index, previous_index, by = 1)
                    index_order <- seq(previous_index, current_index, by = -1)
                    
                    # } else {
                    #   index_order <- seq(previous_index, current_index, by = 1)
                    # }
                  } else {
                    # index_order <- seq(previous_index, current_index, by = 1)
                    index_order <- seq(current_index, previous_index, by = -1)
                  }
                  ### END EDIT
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_snake[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  ############
                  # stepwise_states %>% 
                  #   bind_rows(., implicit_state) -> stepwise_states
                  
                  # 1.26.22 #
                  # Also, change the value in the original dataframe
                  
                  # Remove the last row, insert the row, and add this one on
                  last_row <- stepwise_states[nrow(stepwise_states),]
                  
                  stepwise_states[-nrow(stepwise_states),] %>% 
                    bind_rows(., implicit_state[1,]) %>% 
                    bind_rows(., last_row) -> stepwise_states
                  
                  # stepwise_states <-
                  #   insertRow(stepwise_states, implicit_state[1,], i)
                  
                  ############
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (current_index < previous_index){
                    index_order <- seq(current_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(current_index, previous_index, by = -1)
                  }
                  
                  missing_site <- site_order_notrib_snake[index_order[j+1]]
                  # missing_site <- site_order_notrib_snake[index_order[j]]
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
              }
              
            }
            
            # If it was previously in a tributary
            else if (det_hist[i - 1, 'state'] %in% tributary_mainstem$tributary) {
              # If it was previously in a tributary and this observation isn't
              # in the corresponding mainstem segment, add a line to the detection
              # history
              if (det_hist[i, 'state'] !=  subset(tributary_mainstem, tributary ==
                                                  det_hist[i - 1, 'state'])$mainstem) {
                mainstem_site <-
                  subset(tributary_mainstem, tributary == det_hist[i - 1, 'state'])$mainstem
                
                # INSERT ROW FOR THE MAINSTEM SITE
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = mainstem_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                # stepwise_states %>% 
                #   bind_rows(., implicit_state) -> stepwise_states
                #1.26.22.2:44pm#
                stepwise_states <-
                  insertRow(stepwise_states, implicit_state, i)
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        mainstem_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
                
              }
              
              # If the current state is the correct part of the mainstem, store it
              else if (det_hist[i, 'state'] ==  subset(tributary_mainstem, tributary ==
                                                       det_hist[i - 1, 'state'])$mainstem){
                # Store the tag code
                stepwise_states[i, 'tag_code'] <-
                  det_hist[i, 'tag_code']
                # Store the current state
                stepwise_states[i, 'state'] <-
                  det_hist[i, 'state']
                # Store the time entering this state
                stepwise_states[i, 'date_time'] <-
                  det_hist[i, 'end_time']
                # Store the transition site
                stepwise_states[i, 'pathway'] <-
                  det_hist[i, 'site_class']
              }
              
              ####### THIS CODE COPIED FROM BELOW
              # Zeroth: Jumping between Snake and Columbia, without visiting
              # shared sites in between. We see this with some fish, if they are
              # seen at PRA and then ICH or vice versa, without anything in between
              else if (det_hist[i - 1, 'state'] %in% upper_columbia_sites &
                       det_hist[i, 'state'] %in% snake_sites |
                       det_hist[i - 1, 'state'] %in% snake_sites &
                       det_hist[i, 'state'] %in% upper_columbia_sites) {
                # These inherently have to skip, since they're not seen
                # at a shared site (MCN to ICH/PRA is the only route, and by
                # definition they're not seen there)
                # Get the index of the current site
                current_index <-
                  which(site_order_notrib_columbia_snake %in% det_hist[i, 'state'])
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_columbia_snake %in% det_hist[i - 1, 'state'])
                
                # sequence from current to previous index
                if (current_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, current_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, current_index, by = 1)
                }
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order) - 2)) {
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia_snake[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (current_index < previous_index){
                    index_order <- seq(current_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(current_index, previous_index, by = -1)
                  }
                  
                  missing_site <- site_order_notrib_columbia_snake[index_order[j+1]]
                  # missing_site <- site_order_notrib_columbia_snake[index_order[j]]
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
                
                
              }
              
              
              # First: Shared sites - need two options, depending on where it previously was
              # Shared sites can use either Columbia or Snake sites
              # Condition: If the indices are off by more than one (plus or minus one)
              else if (det_hist[i, 'state'] %in% shared_sites) {
                # If it was previously in the Columbia
                if (det_hist[i - 1, 'state'] %in% c(shared_sites, upper_columbia_sites)) {
                  # Make sure that if it was seen in the adult fishways at a dam, the
                  # detection history contains the state downstream of that dam
                  if (det_hist[i, 'site_class'] %in% c(
                    "BON (adult)",
                    "MCN (adult)",
                    "PRA (adult)",
                    "RIS (adult)",
                    "RRE (adult)",
                    "WEL (adult)",
                    "ICH (adult)",
                    "LGR (adult)"
                  ) &
                  det_hist[i - 1, 'state'] != 
                  site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1]){
                    # (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
                    # which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
                    
                    # Get the missing index
                    missing_index <-
                      which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                    # Get the index of the previous site
                    previous_index <-
                      which(site_order_notrib_snake %in% subset(tributary_mainstem, tributary ==
                                                                  det_hist[i - 1, 'state'])$mainstem)
                    # # Get the missing index
                    # missing_index <-
                    #   which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
                    # # Get the index of the previous site
                    # previous_index <-
                    #   which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
                    
                    # sequence from current to previous index
                    if (missing_index < previous_index) {
                      # index_order <- seq(current_index, previous_index, by = 1)
                      index_order <- seq(previous_index, missing_index, by = -1)
                      
                    } else {
                      index_order <- seq(previous_index, missing_index, by = 1)
                    }
                    
                    # Count the number of sites you need to add and loop through
                    # for (j in 1:(length(index_order) - 2)) {
                    # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                    for (j in 1:(length(index_order) - 1)) {
                      # Add a row for each of them
                      
                      # Insert a new row into stepwise states, with the implicit detection site
                      # Tag code, state, and time (which is NA)
                      missing_site <-
                        site_order_notrib_columbia[index_order[1 + j]]
                      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                   state = missing_site,
                                                   date_time = NA,
                                                   pathway = "implicit")
                      
                      stepwise_states %>% 
                        bind_rows(., implicit_state) -> stepwise_states
                      
                      # Insert a new row into original detection history, with
                      # implicit detection site info
                      
                      # Need to flip the order of sites for these - but it depends on the order of the sites
                      if (missing_index < previous_index){
                        index_order <- seq(missing_index, previous_index, by = 1)
                      }
                      else {
                        index_order <- seq(missing_index, previous_index, by = -1)
                      }
                      
                      # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                      missing_site <- site_order_notrib_columbia[index_order[j]]
                      
                      implicit_detection <- c(det_hist[i, 'tag_code'],
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              "implicit",
                                              missing_site)
                      
                      # Also, change the value in the original dataframe
                      det_hist <-
                        insertRow(det_hist, implicit_detection, i)
                      
                      # Add 1 to the number of added rows
                      added_rows <- added_rows + 1
                    }
                    
                  }
                  
                  else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                           which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) + 1 |
                           which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                           which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) - 1) {
                    # If we see it at two non-consecutive mainstem locations
                    
                    # If the next site skips sites, insert the missing sites
                    # Get the index of the current site
                    current_index <-
                      which(site_order_notrib_columbia %in% det_hist[i, 'state'])
                    # Get the index of the previous site
                    previous_index <-
                      which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
                    
                    # sequence from current to previous index
                    if (current_index < previous_index) {
                      # index_order <- seq(current_index, previous_index, by = 1)
                      index_order <- seq(previous_index, current_index, by = -1)
                      
                    } else {
                      index_order <- seq(previous_index, current_index, by = 1)
                    }
                    
                    # Count the number of sites you need to add and loop through
                    for (j in 1:(length(index_order) - 2)) {
                      # Add a row for each of them
                      
                      # Insert a new row into stepwise states, with the implicit detection site
                      # Tag code, state, and time (which is NA)
                      missing_site <-
                        site_order_notrib_columbia[index_order[1 + j]]
                      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                   state = missing_site,
                                                   date_time = NA,
                                                   pathway = "implicit")
                      
                      stepwise_states %>% 
                        bind_rows(., implicit_state) -> stepwise_states
                      
                      # Insert a new row into original detection history, with
                      # implicit detection site info
                      
                      # Need to flip the order of sites for these - but it depends on the order of the sites
                      if (current_index < previous_index){
                        index_order <- seq(current_index, previous_index, by = 1)
                      }
                      else {
                        index_order <- seq(current_index, previous_index, by = -1)
                      }
                      
                      # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                      missing_site <- site_order_notrib_columbia[index_order[j]]
                      
                      
                      implicit_detection <- c(det_hist[i, 'tag_code'],
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              "implicit",
                                              missing_site)
                      
                      # Also, change the value in the original dataframe
                      det_hist <-
                        insertRow(det_hist, implicit_detection, i)
                      
                      # Add 1 to the number of added rows
                      added_rows <- added_rows + 1
                    }
                  }
                  
                  
                  
                  # Finally, if we're not missing any sites in between, just store it
                  else {
                    # Store the tag code
                    stepwise_states[i, 'tag_code'] <-
                      det_hist[i, 'tag_code']
                    # Store the current state
                    stepwise_states[i, 'state'] <-
                      det_hist[i, 'state']
                    # Store the time entering this state
                    stepwise_states[i, 'date_time'] <-
                      det_hist[i, 'end_time']
                    # Store the transition site
                    stepwise_states[i, 'pathway'] <-
                      det_hist[i, 'site_class']
                  }
                  
                  
                }
                
                # If it was previously in the Snake
                else if (det_hist[i - 1, 'state'] %in% snake_sites) {
                  # Make sure that if it was seen in the adult fishways at a dam, the
                  # detection history contains the state downstream of that dam
                  if (det_hist[i, 'site_class'] %in% c(
                    "BON (adult)",
                    "MCN (adult)",
                    "PRA (adult)",
                    "RIS (adult)",
                    "RRE (adult)",
                    "WEL (adult)",
                    "ICH (adult)",
                    "LGR (adult)"
                  ) &
                  det_hist[i - 1, 'state'] != 
                  site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1]){
                    # (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
                    # which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
                    
                    # Get the missing index
                    missing_index <-
                      which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                    # Get the index of the previous site
                    previous_index <-
                      which(site_order_notrib_snake %in% subset(tributary_mainstem, tributary ==
                                                                  det_hist[i - 1, 'state'])$mainstem)
                    
                    # # Get the missing index
                    # missing_index <-
                    #   which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                    # # Get the index of the previous site
                    # previous_index <-
                    #   which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
                    
                    # sequence from current to previous index
                    if (missing_index < previous_index) {
                      # index_order <- seq(missing_index, previous_index, by = 1)
                      index_order <- seq(previous_index, missing_index, by = -1)
                      
                    } else {
                      index_order <- seq(previous_index, missing_index, by = 1)
                    }
                    
                    # Count the number of sites you need to add and loop through
                    # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                    for (j in 1:(length(index_order) - 1)) {
                      # Add a row for each of them
                      
                      # Insert a new row into stepwise states, with the implicit detection site
                      # Tag code, state, and time (which is NA)
                      missing_site <-
                        site_order_notrib_columbia[index_order[1 + j]]
                      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                   state = missing_site,
                                                   date_time = NA,
                                                   pathway = "implicit")
                      
                      stepwise_states %>% 
                        bind_rows(., implicit_state) -> stepwise_states
                      
                      # Insert a new row into original detection history, with
                      # implicit detection site info
                      
                      # Need to flip the order of sites for these - but it depends on the order of the sites
                      if (missing_index < previous_index){
                        index_order <- seq(missing_index, previous_index, by = 1)
                      }
                      else {
                        index_order <- seq(missing_index, previous_index, by = -1)
                      }
                      
                      # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                      missing_site <- site_order_notrib_columbia[index_order[j]]
                      
                      implicit_detection <- c(det_hist[i, 'tag_code'],
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              "implicit",
                                              missing_site)
                      
                      # Also, change the value in the original dataframe
                      det_hist <-
                        insertRow(det_hist, implicit_detection, i)
                      
                      # Add 1 to the number of added rows
                      added_rows <- added_rows + 1
                    }
                  }
                  
                  else if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                           which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) + 1 |
                           which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                           which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) - 1) {
                    # If we see it at two non-consecutive mainstem locations
                    
                    # If the next site skips sites, insert the missing sites
                    # Get the index of the current site
                    current_index <-
                      which(site_order_notrib_snake %in% det_hist[i, 'state'])
                    # Get the index of the previous site
                    previous_index <-
                      which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
                    
                    # sequence from current to previous index
                    if (current_index < previous_index) {
                      # index_order <- seq(current_index, previous_index, by = 1)
                      index_order <- seq(previous_index, current_index, by = -1)
                      
                    } else {
                      index_order <- seq(previous_index, current_index, by = 1)
                    }
                    
                    # Count the number of sites you need to add and loop through
                    for (j in 1:(length(index_order) - 2)) {
                      # Add a row for each of them
                      
                      # Insert a new row into stepwise states, with the implicit detection site
                      # Tag code, state, and time (which is NA)
                      missing_site <-
                        site_order_notrib_columbia[index_order[1 + j]]
                      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                   state = missing_site,
                                                   date_time = NA,
                                                   pathway = "implicit")
                      
                      stepwise_states %>% 
                        bind_rows(., implicit_state) -> stepwise_states
                      
                      # Insert a new row into original detection history, with
                      # implicit detection site info
                      
                      # Need to flip the order of sites for these - but it depends on the order of the sites
                      if (missing_index < previous_index){
                        index_order <- seq(missing_index, previous_index, by = 1)
                      }
                      else {
                        index_order <- seq(missing_index, previous_index, by = -1)
                      }
                      
                      # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                      missing_site <- site_order_notrib_columbia[index_order[j]]
                      
                      # Insert a new row into original detection history, with
                      # implicit detection site info
                      implicit_detection <- c(det_hist[i, 'tag_code'],
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              NA,
                                              "implicit",
                                              missing_site)
                      
                      # Also, change the value in the original dataframe
                      det_hist <-
                        insertRow(det_hist, implicit_detection, i)
                      
                      # Add 1 to the number of added rows
                      added_rows <- added_rows + 1
                    }
                  }
                  
                  # Finally, if we're not missing any sites in between, just store it
                  else {
                    # Store the tag code
                    stepwise_states[i, 'tag_code'] <-
                      det_hist[i, 'tag_code']
                    # Store the current state
                    stepwise_states[i, 'state'] <-
                      det_hist[i, 'state']
                    # Store the time entering this state
                    stepwise_states[i, 'date_time'] <-
                      det_hist[i, 'end_time']
                    # Store the transition site
                    stepwise_states[i, 'pathway'] <-
                      det_hist[i, 'site_class']
                  }
                  
                }
                
                
              }
              
              # Second: Upper Columbia sites
              else if (det_hist[i, 'state'] %in% upper_columbia_sites) {
                # Make sure that if it was seen in the adult fishways at a dam, the
                # detection history contains the state downstream of that dam
                if (det_hist[i, 'site_class'] %in% c(
                  "BON (adult)",
                  "MCN (adult)",
                  "PRA (adult)",
                  "RIS (adult)",
                  "RRE (adult)",
                  "WEL (adult)",
                  "ICH (adult)",
                  "LGR (adult)"
                ) &
                det_hist[i - 1, 'state'] != 
                site_order_notrib_columbia[which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1]){
                  # (which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1) !=
                  # which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])) {
                  
                  # Get the missing index
                  missing_index <-
                    which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
                  # Get the index of the previous site
                  previous_index <-
                    which(site_order_notrib_columbia %in% subset(tributary_mainstem, tributary ==
                                                                   det_hist[i - 1, 'state'])$mainstem)
                  
                  # sequence from current to previous index
                  if (missing_index < previous_index) {
                    # index_order <- seq(missing_index, previous_index, by = 1)
                    index_order <- seq(previous_index, missing_index, by = -1)
                    
                  } else {
                    index_order <- seq(previous_index, missing_index, by = 1)
                  }
                  
                  # Count the number of sites you need to add and loop through
                  # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                  for (j in 1:(length(index_order) - 1)) {
                    # Add a row for each of them
                    
                    # Insert a new row into stepwise states, with the implicit detection site
                    # Tag code, state, and time (which is NA)
                    missing_site <-
                      site_order_notrib_columbia[index_order[1 + j]]
                    implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                 state = missing_site,
                                                 date_time = NA,
                                                 pathway = "implicit")
                    
                    stepwise_states %>% 
                      bind_rows(., implicit_state) -> stepwise_states
                    
                    # Insert a new row into original detection history, with
                    # implicit detection site info
                    
                    # Need to flip the order of sites for these - but it depends on the order of the sites
                    if (missing_index < previous_index){
                      index_order <- seq(missing_index, previous_index, by = 1)
                    }
                    else {
                      index_order <- seq(missing_index, previous_index, by = -1)
                    }
                    
                    # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                    missing_site <- site_order_notrib_columbia[index_order[j]]
                    
                    implicit_detection <- c(det_hist[i, 'tag_code'],
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            "implicit",
                                            missing_site)
                    
                    # Also, change the value in the original dataframe
                    det_hist <-
                      insertRow(det_hist, implicit_detection, i)
                    
                    # Add 1 to the number of added rows
                    added_rows <- added_rows + 1
                  }
                }
                
                else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                         which(site_order_notrib_columbia %in% det_hist[i -
                                                                        1, 'state']) + 1 |
                         which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                         which(site_order_notrib_columbia %in% det_hist[i -
                                                                        1, 'state']) - 1) {
                  # If we see it at two non-consecutive mainstem locations
                  
                  # If the next site skips sites, insert the missing sites
                  # Get the index of the current site
                  current_index <-
                    which(site_order_notrib_columbia %in% det_hist[i, 'state'])
                  # Get the index of the previous site
                  previous_index <-
                    which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
                  
                  # sequence from current to previous index
                  if (current_index < previous_index) {
                    # index_order <- seq(current_index, previous_index, by = 1)
                    index_order <- seq(previous_index, current_index, by = -1)
                    
                  } else {
                    index_order <- seq(previous_index, current_index, by = 1)
                  }
                  
                  # Count the number of sites you need to add and loop through
                  for (j in 1:(length(index_order) - 2)) {
                    # Add a row for each of them
                    
                    # Insert a new row into stepwise states, with the implicit detection site
                    # Tag code, state, and time (which is NA)
                    missing_site <-
                      site_order_notrib_columbia[index_order[1 + j]]
                    implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                 state = missing_site,
                                                 date_time = NA,
                                                 pathway = "implicit")
                    
                    stepwise_states %>% 
                      bind_rows(., implicit_state) -> stepwise_states
                    
                    # Insert a new row into original detection history, with
                    # implicit detection site info
                    
                    # Need to flip the order of sites for these - but it depends on the order of the sites
                    if (missing_index < previous_index){
                      index_order <- seq(missing_index, previous_index, by = 1)
                    }
                    else {
                      index_order <- seq(missing_index, previous_index, by = -1)
                    }
                    
                    # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                    missing_site <- site_order_notrib_columbia[index_order[j]]
                    
                    implicit_detection <- c(det_hist[i, 'tag_code'],
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            "implicit",
                                            missing_site)
                    
                    # Also, change the value in the original dataframe
                    det_hist <-
                      insertRow(det_hist, implicit_detection, i)
                    
                    # Add 1 to the number of added rows
                    added_rows <- added_rows + 1
                  }
                }
                
                # Finally, if we're not missing any sites in between, just store it
                else {
                  # Store the tag code
                  stepwise_states[i, 'tag_code'] <-
                    det_hist[i, 'tag_code']
                  # Store the current state
                  stepwise_states[i, 'state'] <-
                    det_hist[i, 'state']
                  # Store the time entering this state
                  stepwise_states[i, 'date_time'] <-
                    det_hist[i, 'end_time']
                  # Store the transition site
                  stepwise_states[i, 'pathway'] <-
                    det_hist[i, 'site_class']
                }
                
                
              }
              
              # Third: Snake River sites
              else if (det_hist[i, 'state'] %in% snake_sites) {
                # Make sure that if it was seen in the adult fishways at a dam, the
                # detection history contains the state downstream of that dam
                if (det_hist[i, 'site_class'] %in% c(
                  "BON (adult)",
                  "MCN (adult)",
                  "PRA (adult)",
                  "RIS (adult)",
                  "RRE (adult)",
                  "WEL (adult)",
                  "ICH (adult)",
                  "LGR (adult)"
                ) &
                det_hist[i - 1, 'state'] != 
                site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1]){
                  # (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
                  # which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
                  
                  # Get the missing index
                  missing_index <-
                    which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                  # Get the index of the previous site
                  previous_index <-
                    which(site_order_notrib_snake %in% subset(tributary_mainstem, tributary ==
                                                                det_hist[i - 1, 'state'])$mainstem)
                  
                  # sequence from current to previous index
                  if (missing_index < previous_index) {
                    # index_order <- seq(missing_index, previous_index, by = 1)
                    index_order <- seq(previous_index, missing_index, by = -1)
                    
                  } else {
                    index_order <- seq(previous_index, missing_index, by = 1)
                  }
                  
                  # Count the number of sites you need to add and loop through
                  # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                  for (j in 1:(length(index_order) - 1)) {
                    # Add a row for each of them
                    
                    # Insert a new row into stepwise states, with the implicit detection site
                    # Tag code, state, and time (which is NA)
                    missing_site <-
                      site_order_notrib_columbia[index_order[1 + j]]
                    implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                 state = missing_site,
                                                 date_time = NA,
                                                 pathway = "implicit")
                    
                    stepwise_states %>% 
                      bind_rows(., implicit_state) -> stepwise_states
                    
                    # Insert a new row into original detection history, with
                    # implicit detection site info
                    
                    # Need to flip the order of sites for these - but it depends on the order of the sites
                    if (missing_index < previous_index){
                      index_order <- seq(missing_index, previous_index, by = 1)
                    }
                    else {
                      index_order <- seq(missing_index, previous_index, by = -1)
                    }
                    
                    # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                    missing_site <- site_order_notrib_columbia[index_order[j]]
                    
                    implicit_detection <- c(det_hist[i, 'tag_code'],
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            "implicit",
                                            missing_site)
                    
                    # Also, change the value in the original dataframe
                    det_hist <-
                      insertRow(det_hist, implicit_detection, i)
                    
                    # Add 1 to the number of added rows
                    added_rows <- added_rows + 1
                  }
                }
                
                
                else if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                         which(site_order_notrib_snake %in% det_hist[i -
                                                                     1, 'state']) + 1 |
                         det_hist[i, 'state'] %in% snake_sites &
                         which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                         which(site_order_notrib_snake %in% det_hist[i -
                                                                     1, 'state']) - 1) {
                  # else if (det_hist[i - 1, 'state'] != 
                  #          site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1] |
                  #          det_hist[i - 1, 'state'] != 
                  #          site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) + 1]) {
                  # If we see it at two non-consecutive mainstem locations
                  
                  # If the next site skips sites, insert the missing sites
                  # Get the index of the current site
                  current_index <-
                    which(site_order_notrib_snake %in% det_hist[i, 'state'])
                  # Get the index of the previous site
                  previous_index <-
                    which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
                  
                  # sequence from current to previous index
                  if (current_index < previous_index) {
                    # index_order <- seq(current_index, previous_index, by = 1)
                    index_order <- seq(previous_index, current_index, by = -1)
                    
                  } else {
                    index_order <- seq(previous_index, current_index, by = 1)
                  }
                  
                  # Count the number of sites you need to add and loop through
                  for (j in 1:(length(index_order) - 2)) {
                    # Add a row for each of them
                    
                    # Insert a new row into stepwise states, with the implicit detection site
                    # Tag code, state, and time (which is NA)
                    missing_site <-
                      site_order_notrib_columbia[index_order[1 + j]]
                    implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                                 state = missing_site,
                                                 date_time = NA,
                                                 pathway = "implicit")
                    
                    stepwise_states %>% 
                      bind_rows(., implicit_state) -> stepwise_states
                    
                    # Insert a new row into original detection history, with
                    # implicit detection site info
                    
                    # Need to flip the order of sites for these - but it depends on the order of the sites
                    if (missing_index < previous_index){
                      index_order <- seq(missing_index, previous_index, by = 1)
                    }
                    else {
                      index_order <- seq(missing_index, previous_index, by = -1)
                    }
                    
                    # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                    missing_site <- site_order_notrib_columbia[index_order[j]]
                    
                    implicit_detection <- c(det_hist[i, 'tag_code'],
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            NA,
                                            "implicit",
                                            missing_site)
                    
                    # Also, change the value in the original dataframe
                    det_hist <-
                      insertRow(det_hist, implicit_detection, i)
                    
                    # Add 1 to the number of added rows
                    added_rows <- added_rows + 1
                  }
                }
                
                # Finally, if we're not missing any sites in between, just store it
                else {
                  # Store the tag code
                  stepwise_states[i, 'tag_code'] <-
                    det_hist[i, 'tag_code']
                  # Store the current state
                  stepwise_states[i, 'state'] <-
                    det_hist[i, 'state']
                  # Store the time entering this state
                  stepwise_states[i, 'date_time'] <-
                    det_hist[i, 'end_time']
                  # Store the transition site
                  stepwise_states[i, 'pathway'] <-
                    det_hist[i, 'site_class']
                }
                
                
              }
              
              
              
              
              
              #######
              
              
              
              else {
                # Store the tag code
                stepwise_states[i, 'tag_code'] <-
                  det_hist[i, 'tag_code']
                # Store the current state
                stepwise_states[i, 'state'] <-
                  det_hist[i, 'state']
                # Store the time entering this state
                stepwise_states[i, 'date_time'] <-
                  det_hist[i, 'end_time']
                # Store the transition site
                stepwise_states[i, 'pathway'] <-
                  det_hist[i, 'site_class']
              }
            }
            
          }
          
          # Second: Upper Columbia sites
          else if (det_hist[i, 'state'] %in% upper_columbia_sites) {
            # Make sure that if it was seen in the adult fishways at a dam, the
            # detection history contains the state downstream of that dam
            
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            (which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1) !=
            which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
              
              ### EDIT 2022-07-16: moving this inside for loop
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
              for (j in 1:(length(index_order) - 1)) {
                ### EDIT 2022-07-16
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(missing_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                # end edit
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                     which(site_order_notrib_columbia %in% det_hist[i -
                                                                    1, 'state']) + 1 |
                     which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                     which(site_order_notrib_columbia %in% det_hist[i -
                                                                    1, 'state']) - 1) {
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
              
              # sequence from current to previous index
              if (current_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, current_index, by = -1)
                
              } else {
                # index_order <- seq(previous_index, current_index, by = 1)
                index_order <- seq(current_index, previous_index, by = -1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                ############
                # stepwise_states %>% 
                #   bind_rows(., implicit_state) -> stepwise_states
                
                # Also, change the value in the original dataframe
                stepwise_states <-
                  insertRow(stepwise_states, implicit_state[1,], i)
                
                ############
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            
          }
          
          # Third: Snake River sites
          else if (det_hist[i, 'state'] %in% snake_sites) {
            
            # Make sure that if it was seen in the adult fishways at a dam, the
            # detection history contains the state downstream of that dam
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
            which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
              
              ### EDIT 2022-07-16: moving this inside for loop
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
              for (j in 1:(length(index_order) - 1)) {
                ### EDIT 2022-07-16
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(missing_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                # end edit
                
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            
            else if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                     which(site_order_notrib_snake %in% det_hist[i -
                                                                 1, 'state']) + 1 |
                     which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                     which(site_order_notrib_snake %in% det_hist[i -
                                                                 1, 'state']) - 1) {
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_snake %in% det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
              
              # sequence from current to previous index
              if (current_index < previous_index) {
                index_order <- seq(current_index, previous_index, by = 1)
                # index_order <- seq(previous_index, current_index, by = -1)
                
              # } else {
              #   # index_order <- seq(previous_index, current_index, by = 1)
              #   index_order <- seq(previous_index, current_index, by = -1)
              # }
              } else {
                # index_order <- seq(previous_index, current_index, by = 1)
                index_order <- seq(current_index, previous_index, by = -1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                ############
                # stepwise_states %>% 
                #   bind_rows(., implicit_state) -> stepwise_states
                
                # Also, change the value in the original dataframe
                stepwise_states <-
                  insertRow(stepwise_states, implicit_state[1,], i)
                
                ############
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
          }
          
          ################################################################################################
          
        }          
        
        
        
          

          
          
        
        # If it was in the right part of the mainstem previously, store it like normal
        else {
          # Store the tag code
          stepwise_states[i, 'tag_code'] <-
            det_hist[i, 'tag_code']
          # Store the current state
          stepwise_states[i, 'state'] <-
            det_hist[i, 'state']
          # Store the time entering this state
          stepwise_states[i, 'date_time'] <-
            det_hist[i, 'end_time']
          # Store the transition site
          stepwise_states[i, 'pathway'] <-
            det_hist[i, 'site_class']
        }
        
        
      }
      ### MAINSTEM SITES
      else{
        
        if (det_hist[i - 1, 'state'] %in% tributary_mainstem$tributary) {
          # If it was previously in a tributary and this observation isn't
          # in the corresponding mainstem segment, add a line to the detection
          # history
          if (det_hist[i, 'state'] !=  subset(tributary_mainstem, tributary ==
                                              det_hist[i - 1, 'state'])$mainstem) {
            mainstem_site <-
              subset(tributary_mainstem, tributary == det_hist[i - 1, 'state'])$mainstem
            
            # INSERT ROW FOR THE MAINSTEM SITE
            # Insert a new row into stepwise states, with the implicit detection site
            # Tag code, state, and time (which is NA)
            implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                         state = mainstem_site,
                                         date_time = NA,
                                         pathway = "implicit")
            
            stepwise_states %>% 
              bind_rows(., implicit_state) -> stepwise_states
            
            # Insert a new row into original detection history, with
            # implicit detection site info
            implicit_detection <- c(det_hist[i, 'tag_code'],
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    "implicit",
                                    mainstem_site)
            
            # Also, change the value in the original dataframe
            det_hist <-
              insertRow(det_hist, implicit_detection, i)
            
            # Add 1 to the number of added rows
            added_rows <- added_rows + 1
            
          }
          
          ### EDIT 2022-07-19 ###
          # We have an issue where an individual was seen in the JDR then at the BON (adult) ladders,
          # and the implicit site visit code isn't adding that it must have been in mouth to BON
          # because the order of sites makes sense right now. So we need to first check
          # for detections in adult ladders
          else if (det_hist[i, 'state'] %in% snake_sites) {
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            det_hist[i - 1, 'state'] != 
            site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1]){
              # (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
              # which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% subset(tributary_mainstem, tributary ==
                                                            det_hist[i - 1, 'state'])$mainstem)
              # # Get the missing index
              # missing_index <-
              #   which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
              # # Get the index of the previous site
              # previous_index <-
              #   which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
              
              # MOVE THIS INSIDE FOR LOOP
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              # for (j in 1:(length(index_order) - 2)) {
              # So I don't think we want any minus here, because we're missing both the missing and previous sites
              for (j in 1:(length(index_order))) {
                if (missing_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                # Add a row for each of them
                
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[j]] # make this edit too - I think we can just use j, not j + 1 in this case
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
            
          }
          
          else if (det_hist[i, 'state'] %in% c(upper_columbia_sites, shared_sites)) {
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            det_hist[i - 1, 'state'] != 
            site_order_notrib_columbia[which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1]){
              # (which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1) !=
              # which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% subset(tributary_mainstem, tributary ==
                                                               det_hist[i - 1, 'state'])$mainstem)
              # # Get the missing index
              # missing_index <-
              #   which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
              # # Get the index of the previous site
              # previous_index <-
              #   which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
              
              # MOVE THIS INSIDE FOR LOOP
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              # for (j in 1:(length(index_order) - 2)) {
              # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
              for (j in 1:(length(index_order))) {
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[j]] # make this edit too - I think we can just use j, not j + 1 in this case
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
            
          }
          
          
          ### END 2022-07-19 edit ###
          
          # If the current state is the correct part of the mainstem (and didn't get there via an indirect route), store it
          else if (det_hist[i, 'state'] ==  subset(tributary_mainstem, tributary ==
                                                   det_hist[i - 1, 'state'])$mainstem){
            # Store the tag code
            stepwise_states[i, 'tag_code'] <-
              det_hist[i, 'tag_code']
            # Store the current state
            stepwise_states[i, 'state'] <-
              det_hist[i, 'state']
            # Store the time entering this state
            stepwise_states[i, 'date_time'] <-
              det_hist[i, 'end_time']
            # Store the transition site
            stepwise_states[i, 'pathway'] <-
              det_hist[i, 'site_class']
          }
          
          ####### THIS CODE COPIED FROM BELOW
          # Zeroth: Jumping between Snake and Columbia, without visiting
          # shared sites in between. We see this with some fish, if they are
          # seen at PRA and then ICH or vice versa, without anything in between
          else if (det_hist[i - 1, 'state'] %in% upper_columbia_sites &
                   det_hist[i, 'state'] %in% snake_sites |
                   det_hist[i - 1, 'state'] %in% snake_sites &
                   det_hist[i, 'state'] %in% upper_columbia_sites) {
            # These inherently have to skip, since they're not seen
            # at a shared site (MCN to ICH/PRA is the only route, and by
            # definition they're not seen there)
            # Get the index of the current site
            current_index <-
              which(site_order_notrib_columbia_snake %in% det_hist[i, 'state'])
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_columbia_snake %in% det_hist[i - 1, 'state'])
            
            # sequence from current to previous index
            if (current_index < previous_index) {
              # index_order <- seq(current_index, previous_index, by = 1)
              index_order <- seq(previous_index, current_index, by = -1)
              
            } else {
              index_order <- seq(previous_index, current_index, by = 1)
            }
            
            # Count the number of sites you need to add and loop through
            for (j in 1:(length(index_order) - 2)) {
              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_columbia_snake[index_order[1 + j]]
              implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                           state = missing_site,
                                           date_time = NA,
                                           pathway = "implicit")
              
              stepwise_states %>% 
                bind_rows(., implicit_state) -> stepwise_states
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              
              # Need to flip the order of sites for these - but it depends on the order of the sites
              if (current_index < previous_index){
                index_order <- seq(current_index, previous_index, by = 1)
              }
              else {
                index_order <- seq(current_index, previous_index, by = -1)
              }
              
              missing_site <- site_order_notrib_columbia_snake[index_order[j+1]]
              # missing_site <- site_order_notrib_columbia_snake[index_order[j]]
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              implicit_detection <- c(det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              det_hist <-
                insertRow(det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
            
            
          }
          
          
          # First: Shared sites - need two options, depending on where it previously was
          # Shared sites can use either Columbia or Snake sites
          # Condition: If the indices are off by more than one (plus or minus one)
          else if (det_hist[i, 'state'] %in% shared_sites) {
            # If it was previously in the Columbia
            if (det_hist[i - 1, 'state'] %in% c(shared_sites, upper_columbia_sites)) {
              # Make sure that if it was seen in the adult fishways at a dam, the
              # detection history contains the state downstream of that dam
              if (det_hist[i, 'site_class'] %in% c(
                "BON (adult)",
                "MCN (adult)",
                "PRA (adult)",
                "RIS (adult)",
                "RRE (adult)",
                "WEL (adult)",
                "ICH (adult)",
                "LGR (adult)"
              ) &
              det_hist[i - 1, 'state'] != 
              site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1]){
                # (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
                # which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
                
                # Get the missing index
                missing_index <-
                  which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_snake %in% subset(tributary_mainstem, tributary ==
                                                              det_hist[i - 1, 'state'])$mainstem)
                # # Get the missing index
                # missing_index <-
                #   which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
                # # Get the index of the previous site
                # previous_index <-
                #   which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
                
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                
                # Count the number of sites you need to add and loop through
                # for (j in 1:(length(index_order) - 2)) {
                # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                for (j in 1:(length(index_order) - 1)) {
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (missing_index < previous_index){
                    index_order <- seq(missing_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(missing_index, previous_index, by = -1)
                  }
                  
                  # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                  missing_site <- site_order_notrib_columbia[index_order[j]]
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
                
              }
              
              else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                       which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) + 1 |
                       which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                       which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) - 1) {
                # If we see it at two non-consecutive mainstem locations
                
                # If the next site skips sites, insert the missing sites
                # Get the index of the current site
                current_index <-
                  which(site_order_notrib_columbia %in% det_hist[i, 'state'])
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
                
                # sequence from current to previous index
                if (current_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, current_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, current_index, by = 1)
                }
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order) - 2)) {
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (current_index < previous_index){
                    index_order <- seq(current_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(current_index, previous_index, by = -1)
                  }
                  
                  # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                  missing_site <- site_order_notrib_columbia[index_order[j]]
                  
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
              }
              
              
              
              # Finally, if we're not missing any sites in between, just store it
              else {
                # Store the tag code
                stepwise_states[i, 'tag_code'] <-
                  det_hist[i, 'tag_code']
                # Store the current state
                stepwise_states[i, 'state'] <-
                  det_hist[i, 'state']
                # Store the time entering this state
                stepwise_states[i, 'date_time'] <-
                  det_hist[i, 'end_time']
                # Store the transition site
                stepwise_states[i, 'pathway'] <-
                  det_hist[i, 'site_class']
              }
              
              
            }
            
            # If it was previously in the Snake
            else if (det_hist[i - 1, 'state'] %in% snake_sites) {
              # Make sure that if it was seen in the adult fishways at a dam, the
              # detection history contains the state downstream of that dam
              if (det_hist[i, 'site_class'] %in% c(
                "BON (adult)",
                "MCN (adult)",
                "PRA (adult)",
                "RIS (adult)",
                "RRE (adult)",
                "WEL (adult)",
                "ICH (adult)",
                "LGR (adult)"
              ) &
              det_hist[i - 1, 'state'] != 
              site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1]){
                # (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
                # which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
                
                # Get the missing index
                missing_index <-
                  which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_snake %in% subset(tributary_mainstem, tributary ==
                                                              det_hist[i - 1, 'state'])$mainstem)
                
                # # Get the missing index
                # missing_index <-
                #   which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
                # # Get the index of the previous site
                # previous_index <-
                #   which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
                
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(missing_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                
                # Count the number of sites you need to add and loop through
                # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
                for (j in 1:(length(index_order) - 1)) {
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (missing_index < previous_index){
                    index_order <- seq(missing_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(missing_index, previous_index, by = -1)
                  }
                  
                  # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                  missing_site <- site_order_notrib_columbia[index_order[j]]
                  
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
              }
              
              else if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                       which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) + 1 |
                       which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                       which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) - 1) {
                # If we see it at two non-consecutive mainstem locations
                
                # If the next site skips sites, insert the missing sites
                # Get the index of the current site
                current_index <-
                  which(site_order_notrib_snake %in% det_hist[i, 'state'])
                # Get the index of the previous site
                previous_index <-
                  which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
                
                # sequence from current to previous index
                if (current_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, current_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, current_index, by = 1)
                }
                
                # Count the number of sites you need to add and loop through
                for (j in 1:(length(index_order) - 2)) {
                  # Add a row for each of them
                  
                  # Insert a new row into stepwise states, with the implicit detection site
                  # Tag code, state, and time (which is NA)
                  missing_site <-
                    site_order_notrib_columbia[index_order[1 + j]]
                  implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                               state = missing_site,
                                               date_time = NA,
                                               pathway = "implicit")
                  
                  stepwise_states %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  
                  # Need to flip the order of sites for these - but it depends on the order of the sites
                  if (missing_index < previous_index){
                    index_order <- seq(missing_index, previous_index, by = 1)
                  }
                  else {
                    index_order <- seq(missing_index, previous_index, by = -1)
                  }
                  
                  # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                  missing_site <- site_order_notrib_columbia[index_order[j]]
                  
                  # Insert a new row into original detection history, with
                  # implicit detection site info
                  implicit_detection <- c(det_hist[i, 'tag_code'],
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          NA,
                                          "implicit",
                                          missing_site)
                  
                  # Also, change the value in the original dataframe
                  det_hist <-
                    insertRow(det_hist, implicit_detection, i)
                  
                  # Add 1 to the number of added rows
                  added_rows <- added_rows + 1
                }
              }
              
              # Finally, if we're not missing any sites in between, just store it
              else {
                # Store the tag code
                stepwise_states[i, 'tag_code'] <-
                  det_hist[i, 'tag_code']
                # Store the current state
                stepwise_states[i, 'state'] <-
                  det_hist[i, 'state']
                # Store the time entering this state
                stepwise_states[i, 'date_time'] <-
                  det_hist[i, 'end_time']
                # Store the transition site
                stepwise_states[i, 'pathway'] <-
                  det_hist[i, 'site_class']
              }
              
            }
            
            
          }
          
          # Second: Upper Columbia sites
          else if (det_hist[i, 'state'] %in% upper_columbia_sites) {
            # Make sure that if it was seen in the adult fishways at a dam, the
            # detection history contains the state downstream of that dam
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            det_hist[i - 1, 'state'] != 
            site_order_notrib_columbia[which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1]){
              # (which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1) !=
              # which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% subset(tributary_mainstem, tributary ==
                                                               det_hist[i - 1, 'state'])$mainstem)
              
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
              for (j in 1:(length(index_order) - 1)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                     which(site_order_notrib_columbia %in% det_hist[i -
                                                                    1, 'state']) + 1 |
                     which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                     which(site_order_notrib_columbia %in% det_hist[i -
                                                                    1, 'state']) - 1) {
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
              
              # sequence from current to previous index
              if (current_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, current_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            # Finally, if we're not missing any sites in between, just store it
            else {
              # Store the tag code
              stepwise_states[i, 'tag_code'] <-
                det_hist[i, 'tag_code']
              # Store the current state
              stepwise_states[i, 'state'] <-
                det_hist[i, 'state']
              # Store the time entering this state
              stepwise_states[i, 'date_time'] <-
                det_hist[i, 'end_time']
              # Store the transition site
              stepwise_states[i, 'pathway'] <-
                det_hist[i, 'site_class']
            }
            
            
          }
          
          # Third: Snake River sites
          else if (det_hist[i, 'state'] %in% snake_sites) {
            # Make sure that if it was seen in the adult fishways at a dam, the
            # detection history contains the state downstream of that dam
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            det_hist[i - 1, 'state'] != 
            site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1]){
              # (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
              # which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% subset(tributary_mainstem, tributary ==
                                                            det_hist[i - 1, 'state'])$mainstem)
              
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
              for (j in 1:(length(index_order) - 1)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            
            else if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                     which(site_order_notrib_snake %in% det_hist[i -
                                                                 1, 'state']) + 1 |
                     det_hist[i, 'state'] %in% snake_sites &
                     which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                     which(site_order_notrib_snake %in% det_hist[i -
                                                                 1, 'state']) - 1) {
              # else if (det_hist[i - 1, 'state'] != 
              #          site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1] |
              #          det_hist[i - 1, 'state'] != 
              #          site_order_notrib_snake[which(site_order_notrib_snake %in% det_hist[i, 'state']) + 1]) {
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_snake %in% det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
              
              # sequence from current to previous index
              if (current_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, current_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            # Finally, if we're not missing any sites in between, just store it
            else {
              # Store the tag code
              stepwise_states[i, 'tag_code'] <-
                det_hist[i, 'tag_code']
              # Store the current state
              stepwise_states[i, 'state'] <-
                det_hist[i, 'state']
              # Store the time entering this state
              stepwise_states[i, 'date_time'] <-
                det_hist[i, 'end_time']
              # Store the transition site
              stepwise_states[i, 'pathway'] <-
                det_hist[i, 'site_class']
            }
            
            
          }
          
          
          
          
          
          #######
          
          
          
          else {
            # Store the tag code
            stepwise_states[i, 'tag_code'] <-
              det_hist[i, 'tag_code']
            # Store the current state
            stepwise_states[i, 'state'] <-
              det_hist[i, 'state']
            # Store the time entering this state
            stepwise_states[i, 'date_time'] <-
              det_hist[i, 'end_time']
            # Store the transition site
            stepwise_states[i, 'pathway'] <-
              det_hist[i, 'site_class']
          }
        }
        
        
        # Zeroth: Jumping between Snake and Columbia, without visiting
        # shared sites in between. We see this with some fish, if they are
        # seen at PRA and then ICH or vice versa, without anything in between
        else if (det_hist[i - 1, 'state'] %in% upper_columbia_sites &
                 det_hist[i, 'state'] %in% snake_sites |
                 det_hist[i - 1, 'state'] %in% snake_sites &
                 det_hist[i, 'state'] %in% upper_columbia_sites) {
          # These inherently have to skip, since they're not seen
          # at a shared site (MCN to ICH/PRA is the only route, and by
          # definition they're not seen there)
          # Get the index of the current site
          current_index <-
            which(site_order_notrib_columbia_snake %in% det_hist[i, 'state'])
          # Get the index of the previous site
          previous_index <-
            which(site_order_notrib_columbia_snake %in% det_hist[i - 1, 'state'])
          
          # EDIT 2022-06-10: Move a copy of this ifelse inside the next for loop
          # sequence from current to previous index
          if (current_index < previous_index) {
            # index_order <- seq(current_index, previous_index, by = 1)
            index_order <- seq(previous_index, current_index, by = -1)
            
          } else {
            index_order <- seq(previous_index, current_index, by = 1)
          }
          
          # Count the number of sites you need to add and loop through
          for (j in 1:(length(index_order) - 2)) {
            # sequence from current to previous index
            if (current_index < previous_index) {
              # index_order <- seq(current_index, previous_index, by = 1)
              index_order <- seq(previous_index, current_index, by = -1)
              
            } else {
              index_order <- seq(previous_index, current_index, by = 1)
            }
            # Add a row for each of them
            
            # Insert a new row into stepwise states, with the implicit detection site
            # Tag code, state, and time (which is NA)
            missing_site <-
              site_order_notrib_columbia_snake[index_order[1 + j]]
            implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                         state = missing_site,
                                         date_time = NA,
                                         pathway = "implicit")
            
            stepwise_states %>% 
              bind_rows(., implicit_state) -> stepwise_states
            
            # Insert a new row into original detection history, with
            # implicit detection site info
            
            # Need to flip the order of sites for these - but it depends on the order of the sites
            if (current_index < previous_index){
              index_order <- seq(current_index, previous_index, by = 1)
            }
            else {
              index_order <- seq(current_index, previous_index, by = -1)
            }
            
            missing_site <- site_order_notrib_columbia_snake[index_order[j+1]]
            # missing_site <- site_order_notrib_columbia_snake[index_order[j]]
            
            # Insert a new row into original detection history, with
            # implicit detection site info
            implicit_detection <- c(det_hist[i, 'tag_code'],
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    "implicit",
                                    missing_site)
            
            # Also, change the value in the original dataframe
            det_hist <-
              insertRow(det_hist, implicit_detection, i)
            
            # Add 1 to the number of added rows
            added_rows <- added_rows + 1
          }
          
          
        }
        
        
        # First: Shared sites - need two options, depending on where it previously was
        # Shared sites can use either Columbia or Snake sites
        # Condition: If the indices are off by more than one (plus or minus one)
        else if (det_hist[i, 'state'] %in% shared_sites) {
          # If it was previously in the Columbia
          if (det_hist[i - 1, 'state'] %in% c(shared_sites, upper_columbia_sites)) {
            # Make sure that if it was seen in the adult fishways at a dam, the
            # detection history contains the state downstream of that dam
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            (which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1) !=
            which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
              
              # EDIT: 2022-06-10
              # I'm moving acopy of the next 7 lines inside the for loop below, because otherwise it keeps getting overwritten
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }

              
              # Count the number of sites you need to add and loop through
              # for (j in 1:(length(index_order) - 2)) {
              # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
              for (j in 1:(length(index_order) - 1)) {
                # sequence from current to previous index
                if (missing_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_columbia[index_order[j+1]]
                missing_site <- site_order_notrib_columbia[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
              
            }
            
            else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                     which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) + 1 |
                     which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                     which(site_order_notrib_columbia %in% det_hist[i - 1, 'state']) - 1) {
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
              
              ### EDIT 2022-07-16: moving these lines inside for loop
              # sequence from current to previous index
              if (current_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, current_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                ### EDIT 2022-07-16
                # sequence from current to previous index
                if (current_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, current_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, current_index, by = 1)
                }
                # end edit
                
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>% 
                  bind_rows(., implicit_state) -> stepwise_states
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (current_index < previous_index){
                  index_order <- seq(current_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(current_index, previous_index, by = -1)
                }
                
                #1.26.22.2:13#
                missing_site <- site_order_notrib_columbia[index_order[j+1]]
                # missing_site <- site_order_notrib_columbia[index_order[j]]
                
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            
            
            
            # If we have already inserted the missing site into the stepwise states, ignore it
            # else if (replace_na(stepwise_states[i, 'state'], "NA") == replace_na(stepwise_states[i-1, 'state'], "NA")){
            #   #1.26.22#
            # }
            
            # Finally, if we're not missing any sites in between, 
            # and it's not already in the stepwise states DF, just store it
            else {
              # Store the tag code
              stepwise_states[i, 'tag_code'] <-
                det_hist[i, 'tag_code']
              # Store the current state
              stepwise_states[i, 'state'] <-
                det_hist[i, 'state']
              # Store the time entering this state
              stepwise_states[i, 'date_time'] <-
                det_hist[i, 'end_time']
              # Store the transition site
              stepwise_states[i, 'pathway'] <-
                det_hist[i, 'site_class']
            }
            
            
          }
          
          # If it was previously in the Snake
          else if (det_hist[i - 1, 'state'] %in% snake_sites) {
            # Make sure that if it was seen in the adult fishways at a dam, the
            # detection history contains the state downstream of that dam
            if (det_hist[i, 'site_class'] %in% c(
              "BON (adult)",
              "MCN (adult)",
              "PRA (adult)",
              "RIS (adult)",
              "RRE (adult)",
              "WEL (adult)",
              "ICH (adult)",
              "LGR (adult)"
            ) &
            (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
            which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
              
              # Get the missing index
              missing_index <-
                which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
              
              # Add a row for each of them
              # EDIT: 2022-06-10
              # I'm moving a copy of the next 7 lines inside the for loop below, because otherwise it keeps getting overwritten
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
              for (j in 1:(length(index_order) - 1)) {
                #1.26.22.6:24#
                # Add a row for each of them
                # sequence from current to previous index
                # Refresh this, since index order is changed by the lines below
                if (missing_index < previous_index) {
                  # index_order <- seq(missing_index, previous_index, by = 1)
                  index_order <- seq(previous_index, missing_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, missing_index, by = 1)
                }
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_snake[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                # 1.22.26.6:52 #
                # # # 
                # Need to append zeros, but only if i > nrow
                if (i > nrow(stepwise_states)){
                  NA_rows <- data.frame(tag_code = rep(as.character(NA), i - (nrow(stepwise_states) + 1)),
                                        state = rep(as.character(NA), i - (nrow(stepwise_states) + 1)),
                                        date_time = rep(as.POSIXct(as.character(NA)), i - (nrow(stepwise_states) + 1)),
                                        pathway = rep(as.character(NA), i - (nrow(stepwise_states) + 1)))
                  
                  stepwise_states %>%
                    # Add rows of NAs for any rows that were the same state
                    bind_rows(., NA_rows) %>% 
                    bind_rows(., implicit_state) -> stepwise_states
                }
                
                else{
                  stepwise_states %>%
                    bind_rows(., implicit_state) -> stepwise_states
                }
                # # # 
                
                
                # stepwise_states %>%
                #   bind_rows(., implicit_state) -> stepwise_states
                # Need to run insertRow instead
                # stepwise_states <-
                #   insertRow(stepwise_states, implicit_state, i)
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (missing_index < previous_index){
                  index_order <- seq(missing_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(missing_index, previous_index, by = -1)
                }
                
                # missing_site <- site_order_notrib_snake[index_order[j+1]]
                missing_site <- site_order_notrib_snake[index_order[j]]
                
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            else if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                     which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) + 1 |
                     which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                     which(site_order_notrib_snake %in% det_hist[i - 1, 'state']) - 1) {
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_snake %in% det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
              
              # EDIT: 2022-07-11
              # I'm moving a copy of the next 7 lines inside the for loop below, because otherwise it keeps getting overwritten
            
              # sequence from current to previous index
              if (current_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, current_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                
                if (current_index < previous_index) {
                  # index_order <- seq(current_index, previous_index, by = 1)
                  index_order <- seq(previous_index, current_index, by = -1)
                  
                } else {
                  index_order <- seq(previous_index, current_index, by = 1)
                }
                
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_snake[index_order[1 + j]]
                implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                             state = missing_site,
                                             date_time = NA,
                                             pathway = "implicit")
                
                stepwise_states %>%
                  bind_rows(., implicit_state) -> stepwise_states
                
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                
                # Need to flip the order of sites for these - but it depends on the order of the sites
                if (current_index < previous_index){
                  index_order <- seq(current_index, previous_index, by = 1)
                }
                else {
                  index_order <- seq(current_index, previous_index, by = -1)
                }
                
                missing_site <- site_order_notrib_snake[index_order[j+1]]
                # missing_site <- site_order_notrib_snake[index_order[j]]
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                implicit_detection <- c(det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                det_hist <-
                  insertRow(det_hist, implicit_detection, i)
                
                
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            # Finally, if we're not missing any sites in between, just store it
            else {
              # Store the tag code
              stepwise_states[i, 'tag_code'] <-
                det_hist[i, 'tag_code']
              # Store the current state
              stepwise_states[i, 'state'] <-
                det_hist[i, 'state']
              # Store the time entering this state
              stepwise_states[i, 'date_time'] <-
                det_hist[i, 'end_time']
              # Store the transition site
              stepwise_states[i, 'pathway'] <-
                det_hist[i, 'site_class']
            }
            
          }
          
          
        }
        
        # Second: Upper Columbia sites
        else if (det_hist[i, 'state'] %in% upper_columbia_sites) {
          # Make sure that if it was seen in the adult fishways at a dam, the
          # detection history contains the state downstream of that dam
          if (det_hist[i, 'site_class'] %in% c(
            "BON (adult)",
            "MCN (adult)",
            "PRA (adult)",
            "RIS (adult)",
            "RRE (adult)",
            "WEL (adult)",
            "ICH (adult)",
            "LGR (adult)"
          ) &
          (which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1) !=
          which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])) {
            
            # Get the missing index
            missing_index <-
              which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
            
            # EDIT: 2022-06-10
            # I'm moving a copy of the next 7 lines inside the for loop below, because otherwise it keeps getting overwritten
            # sequence from current to previous index
            if (missing_index < previous_index) {
              # index_order <- seq(missing_index, previous_index, by = 1)
              index_order <- seq(previous_index, missing_index, by = -1)
              
            } else {
              index_order <- seq(previous_index, missing_index, by = 1)
            }
            
            # Count the number of sites you need to add and loop through
            # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
            for (j in 1:(length(index_order) - 1)) {
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              
              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_columbia[index_order[1 + j]]
              implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                           state = missing_site,
                                           date_time = NA,
                                           pathway = "implicit")
              
              stepwise_states %>% 
                bind_rows(., implicit_state) -> stepwise_states
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              
              # Need to flip the order of sites for these - but it depends on the order of the sites
              if (missing_index < previous_index){
                index_order <- seq(missing_index, previous_index, by = 1)
              }
              else {
                index_order <- seq(missing_index, previous_index, by = -1)
              }
              
              # missing_site <- site_order_notrib_columbia[index_order[j+1]]
              missing_site <- site_order_notrib_columbia[index_order[j]]
              
              implicit_detection <- c(det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              det_hist <-
                insertRow(det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
          }
          
          else if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                   which(site_order_notrib_columbia %in% det_hist[i -
                                                                  1, 'state']) + 1 |
                   which(site_order_notrib_columbia %in% det_hist[i, 'state']) <
                   which(site_order_notrib_columbia %in% det_hist[i -
                                                                  1, 'state']) - 1) {
            # If we see it at two non-consecutive mainstem locations
            
            # Get the missing index
            if (which(site_order_notrib_columbia %in% det_hist[i, 'state']) >
                which(site_order_notrib_columbia %in% det_hist[i -
                                                               1, 'state']) + 1){
              missing_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state']) - 1
            }
            else {
              missing_index <-
                which(site_order_notrib_columbia %in% det_hist[i, 'state']) + 1
            }
            
            # If the next site skips sites, insert the missing sites
            # Get the index of the current site
            current_index <-
              which(site_order_notrib_columbia %in% det_hist[i, 'state'])
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_columbia %in% det_hist[i - 1, 'state'])
            
            ### EDIT 2022-07-16: moving these lines inside for loop
            # sequence from current to previous index
            if (current_index < previous_index) {
              # index_order <- seq(current_index, previous_index, by = 1)
              index_order <- seq(previous_index, current_index, by = -1)
              
            } else {
              index_order <- seq(previous_index, current_index, by = 1)
            }
            
            # Count the number of sites you need to add and loop through
            for (j in 1:(length(index_order) - 2)) {
              ### EDIT 2022-07-16
              # sequence from current to previous index
              if (current_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, current_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              # end edit
              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_columbia[index_order[1 + j]]
              implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                           state = missing_site,
                                           date_time = NA,
                                           pathway = "implicit")
              
              stepwise_states %>% 
                bind_rows(., implicit_state) -> stepwise_states
              
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              
              # Need to flip the order of sites for these - but it depends on the order of the sites
              if (missing_index < previous_index){
                index_order <- seq(missing_index, previous_index, by = 1)
              }
              else {
                index_order <- seq(missing_index, previous_index, by = -1)
              }
              
              # missing_site <- site_order_notrib_columbia[index_order[j+1]]
              missing_site <- site_order_notrib_columbia[index_order[j]]
              
              implicit_detection <- c(det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              det_hist <-
                insertRow(det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
          }
          
          # Finally, if we're not missing any sites in between, just store it
          else {
            # Store the tag code
            stepwise_states[i, 'tag_code'] <-
              det_hist[i, 'tag_code']
            # Store the current state
            stepwise_states[i, 'state'] <-
              det_hist[i, 'state']
            # Store the time entering this state
            stepwise_states[i, 'date_time'] <-
              det_hist[i, 'end_time']
            # Store the transition site
            stepwise_states[i, 'pathway'] <-
              det_hist[i, 'site_class']
          }
          
          
        }
        
        # Third: Snake River sites
        else if (det_hist[i, 'state'] %in% snake_sites) {
          # Make sure that if it was seen in the adult fishways at a dam, the
          # detection history contains the state downstream of that dam
          if (det_hist[i, 'site_class'] %in% c(
            "BON (adult)",
            "MCN (adult)",
            "PRA (adult)",
            "RIS (adult)",
            "RRE (adult)",
            "WEL (adult)",
            "ICH (adult)",
            "LGR (adult)"
          ) &
          (which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1) !=
          which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])) {
            
            # Get the missing index
            missing_index <-
              which(site_order_notrib_snake %in% det_hist[i, 'state']) - 1
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
            
            # EDIT: 2022-06-10
            # I'm moving acopy of the next 7 lines inside the for loop below, because otherwise it keeps getting overwritten
            # sequence from current to previous index
            if (missing_index < previous_index) {
              # index_order <- seq(missing_index, previous_index, by = 1)
              index_order <- seq(previous_index, missing_index, by = -1)
              
            } else {
              index_order <- seq(previous_index, missing_index, by = 1)
            }
            
            # Count the number of sites you need to add and loop through
            # Intuitively I think it should be -1 here instead of -2 - because we want to get to the missing site
            for (j in 1:(length(index_order) - 1)) {
              # sequence from current to previous index
              if (missing_index < previous_index) {
                # index_order <- seq(missing_index, previous_index, by = 1)
                index_order <- seq(previous_index, missing_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, missing_index, by = 1)
              }
              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_snake[index_order[1 + j]]
              implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                           state = missing_site,
                                           date_time = NA,
                                           pathway = "implicit")
              
              stepwise_states %>% 
                bind_rows(., implicit_state) -> stepwise_states
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              
              # Need to flip the order of sites for these - but it depends on the order of the sites
              if (missing_index < previous_index){
                index_order <- seq(missing_index, previous_index, by = 1)
              }
              else {
                index_order <- seq(missing_index, previous_index, by = -1)
              }
              
              # missing_site <- site_order_notrib_snake[index_order[j+1]]
              missing_site <- site_order_notrib_snake[index_order[j]]
              
              implicit_detection <- c(det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              det_hist <-
                insertRow(det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
          }
          
          
          else if (which(site_order_notrib_snake %in% det_hist[i, 'state']) >
                   which(site_order_notrib_snake %in% det_hist[i -
                                                               1, 'state']) + 1 |
                   det_hist[i, 'state'] %in% snake_sites &
                   which(site_order_notrib_snake %in% det_hist[i, 'state']) <
                   which(site_order_notrib_snake %in% det_hist[i -
                                                               1, 'state']) - 1) {
            # If we see it at two non-consecutive mainstem locations
            
            # If the next site skips sites, insert the missing sites
            # Get the index of the current site
            current_index <-
              which(site_order_notrib_snake %in% det_hist[i, 'state'])
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_snake %in% det_hist[i - 1, 'state'])
            
            ### EDIT 2022-07-16: moving these lines inside next for loop
            # sequence from current to previous index
            if (current_index < previous_index) {
              # index_order <- seq(current_index, previous_index, by = 1)
              index_order <- seq(previous_index, current_index, by = -1)
              
            } else {
              index_order <- seq(previous_index, current_index, by = 1)
            }
            
            # Count the number of sites you need to add and loop through
            for (j in 1:(length(index_order) - 2)) {
              ### EDIT 2022-07-16
              # sequence from current to previous index
              if (current_index < previous_index) {
                # index_order <- seq(current_index, previous_index, by = 1)
                index_order <- seq(previous_index, current_index, by = -1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              # end edit
              
              
              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_snake[index_order[1 + j]]
              
              missing_index <- index_order[1 + j]
              implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                           state = missing_site,
                                           date_time = NA,
                                           pathway = "implicit")
              
              stepwise_states %>% 
                bind_rows(., implicit_state) -> stepwise_states
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              
              # Need to flip the order of sites for these - but it depends on the order of the sites
              if (missing_index < previous_index){
                index_order <- seq(missing_index, previous_index, by = 1)
              }
              else {
                index_order <- seq(missing_index, previous_index, by = -1)
              }
              
              # missing_site <- site_order_notrib_snake[index_order[j+1]]
              missing_site <- site_order_notrib_snake[missing_index]
              
              implicit_detection <- c(det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              det_hist <-
                insertRow(det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
          }
          
          # Finally, if we're not missing any sites in between, just store it
          else {
            # Store the tag code
            stepwise_states[i, 'tag_code'] <-
              det_hist[i, 'tag_code']
            # Store the current state
            stepwise_states[i, 'state'] <-
              det_hist[i, 'state']
            # Store the time entering this state
            stepwise_states[i, 'date_time'] <-
              det_hist[i, 'end_time']
            # Store the transition site
            stepwise_states[i, 'pathway'] <-
              det_hist[i, 'site_class']
          }
          
          
        }
        
        
        # Finally, if we're not missing any sites in between, just store it
        else {
          # Store the tag code
          stepwise_states[i, 'tag_code'] <-
            det_hist[i, 'tag_code']
          # Store the current state
          stepwise_states[i, 'state'] <-
            det_hist[i, 'state']
          # Store the time entering this state
          stepwise_states[i, 'date_time'] <-
            det_hist[i, 'end_time']
          # Store the transition site
          stepwise_states[i, 'pathway'] <-
            det_hist[i, 'site_class']
        }
      }
      
    }
    
    
    
    
    ### REPEAT DETECTION AT THE SAME DAM
    # If it's next seen at the same dam, add a line to indicate it fell back
    # I think we need a statement for each dam, because each is unique
    # EDITS: 06-10-2022
    # So, I'm changing the first line from det_hist[i, 'site_class'] == det_hist[i - 1, 'site_class'] to:
    # det_hist[i, 'state'] == det_hist[i - 1, 'state']
    # This is because sometimes we see fish in the adult ladders after they came
    # from the upstream state, but weren't last in the adult ladders. For example, the inriver arrays.
    # So, by changing this from site class to state, we should make this more general
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "BON (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      
      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, mouth to BON",
                                   date_time = NA,
                                   pathway = "implicit")
      
      # EDIT 06-10-22
      # I'm changing from insertRow to bind_row - in every other instance I think we're using bind_row,
      # and using insertRow here is messing it up. And in each of the following dam lines
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit",
                              "mainstem, mouth to BON")
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "MCN (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, BON to MCN",
                                   date_time = NA,
                                   pathway = "implicit")
      
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit",
                              "mainstem, BON to MCN")
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "PRA (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, MCN to ICH or PRA",
                                   date_time = NA,
                                   pathway = "implicit")
      
      
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(
        det_hist[i, 'tag_code'],
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        "implicit",
        "mainstem, MCN to ICH or PRA"
      )
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "RIS (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, PRA to RIS",
                                   date_time = NA,
                                   pathway = "implicit")
      
      
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit",
                              "mainstem, PRA to RIS")
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "RRE (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, RIS to RRE",
                                   date_time = NA,
                                   pathway = "implicit")
      
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit",
                              "mainstem, RIS to RRE")
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "WEL (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, RRE to WEL",
                                   date_time = NA,
                                   pathway = "implicit")
      
      
      
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit",
                              "mainstem, RRE to WEL")
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "ICH (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, MCN to ICH or PRA",
                                   date_time = NA,
                                   pathway = "implicit")
      
      
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(
        det_hist[i, 'tag_code'],
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        "implicit",
        "mainstem, MCN to ICH or PRA"
      )
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (det_hist[i, 'state'] == det_hist[i - 1, 'state'] &
             det_hist[i, 'site_class'] == "LGR (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)

      implicit_state <- data.frame(tag_code = det_hist[i, 'tag_code'],
                                   state = "mainstem, ICH to LGR",
                                   date_time = NA,
                                   pathway = "implicit")
      
      stepwise_states %>% 
        bind_rows(., implicit_state) -> stepwise_states
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit",
                              "mainstem, ICH to LGR")
      
      # Also, change the value in the original dataframe
      det_hist <- insertRow(det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }

    
    
    # If it's in the same state, don't record it
    else {
      # Nothing!
    }
  }
  # If it's a different fish:
  else {
    # Increase fish counter
    fish_counter <- fish_counter + 1
    print(paste0("Fish #", fish_counter))
    # Start a new entry
    # Store the tag code
    stepwise_states[i, 'tag_code'] <-
      det_hist[i, 'tag_code']
    # Store the current state
    stepwise_states[i, 'state'] <- det_hist[i, 'state']
    # Store the time entering this state
    stepwise_states[i, 'date_time'] <-
      det_hist[i, 'end_time']
    # Store the transition site
    stepwise_states[i, 'pathway'] <-
      det_hist[i, 'site_class']
  }
}


# stepwise_states







# # # Remove all NAs
stepwise_states_noNA <- stepwise_states[!is.na(stepwise_states$tag),]

stepwise_states_noNA  %>% 
  mutate(prev_tag_code = lag(tag_code)) %>% 
  mutate(next_tag_code = lead(tag_code)) %>% 
  mutate(next_state = lead(state)) %>% 
  mutate(next_state = ifelse(next_tag_code != tag_code, "loss", next_state)) -> transitions

# write.csv(transitions, here::here("figures", "bad_transitions_table_2022-07-16.csv"))



# Compare them to original detection histories

# Export this for fitting multistate
write.csv(stepwise_states_noNA, "states_complete.csv")

##### Sort individuals into run years #####
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22")
run_year_start <- seq(ymd("2004-06-01"), ymd("2021-06-01"), by = "years")
run_year_end <- seq(ymd("2005-05-31"), ymd("2022-05-31"), by = "years")

run_year_df <- data.frame(run_year, run_year_start, run_year_end)

det_hist %>%
  group_by(tag_code) %>%
  filter(row_number() == 1) %>%
  mutate(run_year = subset(run_year_df, run_year_start < start_time & run_year_end > start_time)$run_year) %>%
  dplyr::select(tag_code, run_year) -> tag_codes_run_year

stepwise_states_noNA %>%
  left_join(., tag_codes_run_year, by = "tag_code") -> stepwise_states_noNA

# Export this for fitting multistate
write.csv(stepwise_states_noNA, "states_complete.csv")

