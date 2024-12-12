library(sf)
library(dplyr)
library(terra)
library(ggplot2)
sf::sf_use_s2(FALSE)
# Load extracted data
extracted_data <- read.csv("Data/cholera_climate_systematic_review_extract_into_sheets_shps.csv")

# Load raster and convert to polygons
climatic_zones_raster <- terra::rast("Data/climatic_zone_raster/koppen_geiger_0p00833333.tif")
climatic_zones_sf <- sf::st_as_sf(terra::as.polygons(climatic_zones_raster))

# Load shapefile and prepare for join
climatic_shp <- readRDS("Data/shapefiles/location_shapefiles.rds")
climatic_shp <- climatic_shp %>% mutate(lctn_pr = as.character(lctn_pr))

# Estimate the intersections between shapefiles and climatic zone shapefiles
#for (lp_idx in unique(extracted_data)) {
#  
#  shp_tmp <- climatic_shp[climatic_shp$lctn_pr == lp_idx,]$geometry
#  
#  # estimate the intersections
#  intersections_tmp <- sf::st_intersection(climatic_zones_sf, shp_tmp) %>% 
#    mutate(lctn_pr = lp_idx)
#  
#  intersections <- rbind(intersections, intersections_tmp)
#  
#}

#saveRDS(intersections, "Data/shapefiles/intersected_shp.rds")

# Load climatic zone dictionary
climatic_zone_dic <- read.csv("Data/climatic_zone_raster/climatic_zone_dictionary.csv")

# Perform spatial intersection
largest_intersections <- readRDS("Data/shapefiles/intersected_shp.rds") %>% 
  mutate(intersected_area = sf::st_area(geometry)) %>%
  group_by(lctn_pr) %>%
  filter(intersected_area == max(intersected_area, na.rm = TRUE)) %>%
  ungroup() %>% 
  left_join(climatic_zone_dic, by =c("koppen_geiger_0p00833333"="climatic_zone_idx")) 


# Select the relevant columns
result <- largest_intersections %>%
  select(lctn_pr, climatic_zone) %>% 
  as.data.frame() %>% 
  select(!geometry) %>% 
  right_join(extracted_data, by = c("lctn_pr" = "location_period_id")) %>% 
  left_join(climatic_shp, by = "lctn_pr") %>% 
  mutate(climatic_zone = as.factor(climatic_zone)) %>% 
  st_as_sf()


# Load country shapefiles
country_shp <- readRDS("Data/shapefiles/country_shapefiles.rds")

# plot
p_climatic_zone <- ggplot() +
  geom_sf(data=result, aes(fill = climatic_zone),alpha = 0.3) +
  geom_sf(data=country_shp,alpha=0) +
  facet_wrap(~eff_admin, ncol = 2) +
  theme_minimal() 

ggsave("Analysis/Figures/climatic_zone_by_admin_new.png", plot = p_climatic_zone, width = 8, height = 6, dpi = 300)
