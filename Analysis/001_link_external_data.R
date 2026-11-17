# Load relevant R packages ----
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(osmdata)
library(patchwork)
sf::sf_use_s2(FALSE)

# Extracted data ----
extracted_data <- read.csv("Data/extracted_data/cholera_climate_systematic_review_extract_into_sheets_shps.csv")

# Climatic zone ----
## Load raster and convert to polygons
climatic_zones_raster <- terra::rast("Data/climatic_zone_raster/koppen_geiger_0p00833333.tif")
climatic_zones_sf <- sf::st_as_sf(terra::as.polygons(climatic_zones_raster))

## Load shapefile and prepare for join
climatic_shp <- readRDS("Data/shapefiles/location_shapefiles.rds") %>% 
  dplyr::mutate(lctn_pr = as.character(lctn_pr))

## Perform spatial intersection
### Estimate the intersections between shapefiles and climatic zone shapefiles (only need to run one time)
if(!file.exists("Data/shapefiles/intersected_shp.rds")){
  intersections <- NULL
  for (lp_idx in unique(as.character(extracted_data$location_period_id))) {
    
    shp_tmp <- climatic_shp[climatic_shp$lctn_pr == lp_idx,]$geometry
    
    # estimate the intersections
    intersections_tmp <- sf::st_intersection(climatic_zones_sf, shp_tmp) %>%
      dplyr::mutate(lctn_pr = lp_idx)
    
    intersections <- dplyr::bind_rows(intersections, intersections_tmp)
    
  }
  saveRDS(intersections, "Data/shapefiles/intersected_shp.rds")
}

## Load climatic zone dictionary
climatic_zone_dic <- read.csv("Data/climatic_zone_raster/climatic_zone_dictionary.csv")

## Subset the largest intersected areas for each location period id
largest_intersections <- readRDS("Data/shapefiles/intersected_shp.rds") %>% 
  dplyr::mutate(intersected_area = sf::st_area(geometry)) %>%
  dplyr::group_by(lctn_pr) %>%
  dplyr::filter(intersected_area == max(intersected_area, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(climatic_zone_dic, by =c("koppen_geiger_0p00833333"="climatic_zone_idx")) 

## Select the relevant columns
result_climate <- largest_intersections %>%
  dplyr::select(lctn_pr, new_climatic_zone) %>% 
  as.data.frame() %>% 
  dplyr::select(!geometry) %>% 
  dplyr::right_join(extracted_data, by = c("lctn_pr" = "location_period_id")) %>% 
  dplyr::left_join(climatic_shp, by = "lctn_pr") %>% 
  dplyr::mutate(new_climatic_zone = as.factor(new_climatic_zone)) %>% 
  sf::st_as_sf()

## Remove geometry of export_result
export_result <- result_climate %>% 
  sf::st_drop_geometry() 

# Inland and Costcal ----
total_shp<-readRDS("Data/shapefiles/location_shapefiles.rds")

## Historical coastlines data downloaded
coastlines <- sf::st_read("Data/inland_coast/coastlines-split-4326/lines.shp")
#https://osmdata.openstreetmap.de/data/coastlines.html

## Query coastline data
identify_dist_to_water <- function(shp,dist_threshold = 10^5){

  coast <- tryCatch(
    {
      osmdata::opq(bbox = sf::st_bbox(shp)) %>%
      osmdata::add_osm_feature(key = "natural", value = "coastline") %>%
      osmdata::osmdata_sf()
    },
    error = function(e) {
      cat("Primary query failed",e$message,"\n")
      tryCatch(
        {
          sf::st_crop(coastlines,sf::st_bbox(shp))
        }
      )
    }
  )

  if(is.null(coast$osm_lines) & any(names(coast) == "osm_multilines")){
    
    type = "Inland" # not even close to the coast
    
  } else{
    if(!is.null(coast$osm_lines)){
      coast_sf <- coast$osm_lines
    } else {
      coast_sf <- coast
    }
    
    # Set the crs consistent between the coast line and shapefile
    shp <- sf::st_transform(shp %>% sf::st_make_valid(), crs = sf::st_crs(coast_sf))
    
    # Find the cloest coast line point near the shp
    nearest_coastline_segment <- sf::st_nearest_feature(shp, coast_sf)
    
    # Calculate distances to the closest coastline segment
    distances_nearest <- sf::st_distance(shp, coast_sf[nearest_coastline_segment, ])
    
    if (as.numeric(distances_nearest) < dist_threshold) {
      type <- "Coastal"
    } else {
      type <- "Inland"
    }
  }
  return(type)
}

## Estimate the distance to water and decide whether it's inland or costcal
if(!file.exists("Data/inland_coast/dist_to_water_results.rds")){
  dist_to_water_results <- purrr::map(total_shp$geometry, function(geom) {
    single_sf <- sf::st_sf(geometry = sf::st_sfc(geom))
    sf::st_crs(single_sf) <- 4326
    identify_dist_to_water(single_sf)
  })
  saveRDS(dist_to_water_results,"Data/inland_coast/dist_to_water_results.rds")  
} else {
  dist_to_water_results <- readRDS("Data/inland_coast/dist_to_water_results.rds")
}

total_shp$inland_coast <- unlist(dist_to_water_results)

## Perform the left join to add the inland_coast column
export_result <- export_result %>%
  dplyr::left_join(total_shp[, c("lctn_pr", "inland_coast")] %>% sf::st_drop_geometry(), by = "lctn_pr")

# HDI ----
hdi_raster <- raster::raster("Data/hdi_raster/Sub_National_HDI.tif")

## Extract the average Sub_National_HDI for each polygon in total_shp
result_hdi <- exactextractr::exact_extract(hdi_raster, total_shp, fun = "mean", append_cols = TRUE)

## Join the hdi result to exported data
export_result <- export_result %>% 
  dplyr::left_join(result_hdi,by="lctn_pr") %>% 
  dplyr::rename(average_hdi = mean)

# A new spatial scale categorization: national, supranational and subnational ----
extracted_data <-export_result %>% 
  dplyr::mutate(hdi_cat = ifelse(average_hdi < 0.55, "low hdi", 
                                 ifelse(average_hdi < 0.7, "medium hdi", 
                                        ifelse(average_hdi < 0.8, "high hdi", "very high hdi")))) %>% 
  dplyr::mutate(hdi_cat = factor(hdi_cat, levels = c("low hdi", "medium hdi", "high hdi", "very high hdi"))) %>% 
  dplyr::mutate(
    new_admin = case_when(
      eff_admin == "non-adm geography" & stringr::str_detect(eff_loc, "multihospitals|multiprovince|multidistrict|regions|multiregion|Unguja|Pemba neighborhoods\\(shehias\\)|multi health districts|multicommune|Lokoja|CMCH|bengal delta region") & 
        !(eff_country %in% c("IND", "PAK") & eff_loc == "multiregion") ~ "subnational",
      studykey == "Ruiz-Moreno2007" ~ "subnational",
      eff_admin %in% c("adm1", "adm2", "adm3") ~ "subnational",
      eff_admin == "adm0" ~ "national",
      eff_admin == "non-adm geography" & eff_loc == "bengal delta region" ~ "subnational",
      eff_country == "HTI" & eff_admin == "non-adm geography" ~ "subnational",
      eff_country == "NGA" & eff_loc == "multistate" ~ "subnational",
      studykey == "Oloukoi2014" ~ "subnational",
      eff_country == "IND|PAK" & studykey == "Jutla2013" ~ "supranational",
      TRUE ~ "supranational"
    )
  )
#write.csv(extracted_data,"Data/processed_cholera_climate_systematic_review_extract_into_sheets.csv",row.names=F)

# Endemicity ----
## Processed extracted data ----
processed_extracted_data <- read.csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv") %>% 
  filter(new_admin != "supranational") %>% 
  rowwise() %>% 
  mutate(
    year_tl = ifelse(nchar(eff_tl) == 5, as.numeric(lubridate::year(as.Date(as.numeric(eff_tl),origin="1899-12-30"))), as.numeric(substr(eff_tl,1,4))),
    year_tr = ifelse(nchar(eff_tr) == 5, as.numeric(lubridate::year(as.Date(as.numeric(eff_tr),origin="1899-12-30"))), as.numeric(substr(eff_tr,1,4)))
  )

## psql code to pull WHO annual report data for these countries ----
### \copy (select time_left,time_right,suspected_cases,deaths,location_id from observations where (observation_collection_id <=480) and (observation_collection_id>=300) or (observation_collection_id=620)) to '/home/qulu/who_annual_report.csv' (format csv, delimiter ';');

report<-read.csv("Data/endemicity/who_annual_report.csv") %>% 
  tidyr::separate(obs,c('time_left','time_right','suspected_cases','deaths','location_id'),sep = ';',remove=T)
loc<-read.csv("Data/endemicity/country_location_id_name.csv")

report<-merge(report,loc,by.x="location_id",by.y="id",all.x=TRUE) %>% 
  filter(is.na(region)==F) # two OCs are not at the country-level between 300-480

## Create a new column for the midpoint of time
report$time <- (lubridate::ymd(report$time_right) -lubridate::ymd(report$time_left)) / 2 + lubridate::ymd(report$time_left)
report<-report%>% arrange(region)%>%mutate(region = as.factor(region),suspected_cases=as.numeric(suspected_cases))

## Convert year_tl and year_tr to Date format
processed_extracted_data <- processed_extracted_data %>%
  mutate(year_tl = lubridate::ymd(paste0(year_tl, "-01-01")),
         year_tr = lubridate::ymd(paste0(year_tr, "-12-31")))

## Merge region information with study period data
processed_extracted_data <- processed_extracted_data %>%
  left_join(report %>% select(region) %>% distinct(), by = c("eff_country" = "region")) %>% 
  mutate(eff_country = as.factor(eff_country))

report <- report %>% 
  filter(region %in% processed_extracted_data$eff_country) %>% 
  mutate(region = factor(region,levels = levels(processed_extracted_data$eff_country)))

country_positions <- as.numeric(as.factor(report$region))

## Identify countries whose study period TL is earlier than the earliest reported WHO cases
earliest_data <- report %>%
  group_by(region) %>%
  summarise(earliest_year = min(time_left, na.rm = TRUE))

## Merge with processed_extracted_data to compare TL with earliest available data
processed_extracted_data <- processed_extracted_data %>%
  left_join(earliest_data, by = c("eff_country" = "region")) %>%
  mutate(border_color = ifelse(year_tl < earliest_year, "Before the earliest WHO report", "black"))

## Heatmap
p <- ggplot() +
  geom_rect(data=report %>% mutate(suspected_cases = round(suspected_cases,0)),
            aes(y = region, fill = suspected_cases,
                xmin = lubridate::ymd(time_left), xmax = lubridate::ymd(time_right), 
                ymin = as.numeric(as.factor(region)) - 0.5, 
                ymax = as.numeric(as.factor(region)) + 0.5)) +
  scale_fill_viridis_c(option = "D", trans = "log",
                       breaks = c(1, 10, 100, 1000, 10000,100000),  # Adjust breaks
                       labels = c("1", "10", "100", "1000", "10000","100000")) +
  labs(
    x = "Time Interval",
    y = "Country",
    fill ="Annual Cases",
    title = "Heatmap of WHO Reported Cases by Time Interval and Country"
  ) +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "5 year") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
    ) +
  # Add study period highlighting per country
  geom_rect(data = processed_extracted_data %>% distinct(year_tl, year_tr, eff_country, border_color), 
            aes(xmin = year_tl, xmax = year_tr, 
                ymin = as.numeric(as.factor(eff_country)) - 0.5, 
                ymax = as.numeric(as.factor(eff_country)) + 0.5, 
                color = border_color),  # Dynamic border color
            fill = "#FFA500",
            alpha = 0.4) +  
  scale_color_identity(guide = "legend", name = "Study Period", labels = c("After earliest WHO report","Before earliest WHO report")) + 
  geom_hline(yintercept = country_positions, color = "black", linetype = "dotted", alpha = 0.5)
p

#ggplot2::ggsave("Analysis/Figures/Heatmap of WHO annual cases.png", plot = p,width = 10,height=6)

## Function to identify endemicity status of the country based on the years prior to study period ----

#' @name identify_endemicity
#' @export 
#' @param processed_data extracted and processed data
#' @param study_year_col the column name for the year of the study/estimate 
#' @param study_country_col the column name for the country of the study/estimate 
#' @param who_data WHO annual report data
#' @param who_year_col the column name for the year of the WHO annual reports 
#' @param who_country_col the column name for the country of the WHO annual reports
#' @param who_case_col the colname name for the cases of the WHO annual reports
#' @param prior_year_with_sCh the number of years with WHO reported cases before the study year used to define endemicity  
#' @returns processed_data with a new column called endemicity with three levels: Endemic, Non-Endemic, and Unknown

identify_endemicity <- function(
    processed_data = processed_extracted_data[1,] %>%
      mutate(year_tl = lubridate::ymd(year_tl)),
    study_year_col = "year_tl",
    study_country_col = "eff_country",
    who_data = report %>% mutate(year = lubridate::year(time_left)),
    who_year_col = "year",
    who_country_col = "region",
    who_case_col = "suspected_cases",
    prior_year_with_sCh = 3
) {
  
  library(dplyr)
  library(tidyr)
  
  processed_data <- processed_data %>% 
    mutate(!!study_year_col := as.numeric(!!sym(study_year_col)))
  
  processed_data <- processed_data %>%
    rowwise() %>%
    mutate(
      endemicity_status = {
        
        study_year <- get(study_year_col)
        study_country <- get(study_country_col)
        
        prior_years <- seq(study_year - prior_year_with_sCh, study_year - 1)
        
        # Create a full dataset of years and join with WHO data
        full_years_df <- tibble(
          !!who_year_col := prior_years,
          !!who_country_col := study_country
        )
        
        who_filtered <- full_years_df %>%
          left_join(who_data, by = c(who_year_col, who_country_col))
        
        # Check if all years have NA (no reports at all)
        if (all(is.na(who_filtered[[who_case_col]]))) {
          status <- "Unknown"
        } else {
          # Replace NA with 0 if at least one year has a report
          who_filtered[[who_case_col]][is.na(who_filtered[[who_case_col]])] <- 0
          
          reported_years <- sum(who_filtered[[who_case_col]] > 0)
          
          # Determine endemicity status
          status <- ifelse(reported_years >= prior_year_with_sCh, "Endemic", "Non-Endemic")
        }
        status
      }
    ) %>%
    ungroup()
  
  return(processed_data)
}

## proportion of years with reported WHO cases between 2000-2023 ----
library(tidyverse)
library(purrr)
source(here::here("Analysis/utils.R"))
library(gtsummary)
library(lubridate)
library(gt)
library(survey)
library(dplyr)

# Load the dataset from the file path
data <- readr::read_csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv") %>% 
  mutate(
    admin_level = case_when(
      new_admin == "supranational" ~ "regional",
      TRUE ~ new_admin
    )
  ) %>% 
  add_continent("eff_country") %>% 
  mutate(
    exp_class =factor(exp_class, levels = c("precipitation","temperature","flood","drought"))
  ) %>% 
  mutate(across(c("eff_tl", "eff_tr"), ~ char_to_date(.) %>% as.Date())) %>%
  mutate(study_year_tl = lubridate::year(eff_tl),
         study_year_tr = lubridate::year(eff_tr)) 

who_data<-read.csv("Data/endemicity/who_annual_report.csv") %>%
  tidyr::separate(obs,c('time_left','time_right','suspected_cases','deaths','location_id','primary'),sep = ';',remove=T) %>%
  filter(primary == 't') %>% 
  merge(.,read.csv("Data/endemicity/country_location_id_name.csv"),by.x="location_id",by.y="id",all.x=TRUE) %>%
  filter(is.na(region)==F) %>%
  arrange(region) %>%
  mutate(country = as.factor(region),suspected_cases=as.numeric(suspected_cases)) %>% 
  filter(time_right >= as.Date("2000-01-01") & time_right <= as.Date("2023-12-31")) %>% 
  mutate(year = lubridate::year(time_right)) %>% 
  select(country,suspected_cases,year)

endemicity_data <- identify_endemicity(
  processed_data = data,
  who_data = who_data,
  group_col = 'exp_class'
) %>%
  mutate(endemicity_group = case_when(
    percent_year_with_cases <= 50 ~ "0-50%",
    percent_year_with_cases > 50 & percent_year_with_cases <= 95 ~ "50-95%",
    percent_year_with_cases > 95 ~ ">95%"
  ))

data_now <- readr::read_csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv")

data_expanded <- data_now %>% 
  mutate(row_id = row_number()) %>% 
  separate_rows(eff_country, sep = "\\|") %>%
  rename(country = eff_country)

data_joined <- data_expanded %>%
  left_join(endemicity_data %>%
              select(country, exp_class, percent_year_with_cases),
            by = c("country", "exp_class"))

data_summarized <- data_joined %>%
  group_by(row_id) %>%
  summarise(
    eff_country = str_c(unique(country), collapse = "|"),
    exp_class = first(exp_class),
    mean_percent_year_with_cases = mean(percent_year_with_cases, na.rm = TRUE)
  )

data_final <- data_now %>%
  mutate(row_id = row_number()) %>%
  left_join(data_summarized %>% select(row_id,mean_percent_year_with_cases), by = "row_id") %>%
  select(-row_id)
#write.csv(data_final,"Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv",row.names = F)


# Urbanicity ----
library(OutbreakExtractR)
sf::sf_use_s2(FALSE)
pre_data <- read.csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv") %>% 
  mutate(across(c("eff_tl", "eff_tr"), ~ char_to_date(.) %>% as.Date())) %>%
  mutate(study_year_tl = lubridate::year(eff_tl),
         study_year_tr = lubridate::year(eff_tr)) %>% 
  filter(new_admin != 'supranational') %>%  # only add to subantional and national
  filter_final_analysis_dataset(
    ., criteria = list(expr(study_year_tl > 1961))
  ) %>% 
  distinct(
    lctn_pr,
    eff_country,
    study_year_tl
  ) %>% 
  mutate(
    pop = NA,
    area_per_1km2 = NA,
    rural_urban = NA
  )

total_shp <- readRDS("Data/shapefiles/location_shapefiles.rds")
for (row_idx in 1:nrow(pre_data)) {
  lp_id = pre_data$lctn_pr[row_idx]
  year_tmp = ifelse(pre_data$study_year_tl[row_idx] < 2000, 2000, pre_data$study_year_tl[row_idx])
  eff_country <- pre_data$eff_country[row_idx]
  
  country_shp <- geodata::gadm(country = eff_country, level = 0, path = tempdir()) %>% 
    sf::st_as_sf() %>% 
    filter(GID_0 == eff_country)
  
  pop <- get_pop(
    shp = total_shp[total_shp$lctn_pr==lp_id,],
    year = year_tmp,
    country = eff_country,
    country_shp = country_shp,
    pop_raster_path = "Data/world_pop"
  ) 
  
  if(pre_data$study_year_tl[row_idx] < 2000){
    pre_2000_pop <- WPP2022 %>% filter(ISO3_code == eff_country & Time == pre_data$study_year_tl[row_idx]) %>% pull(PopTotal)
    pop_2000 <- WPP2022 %>% filter(ISO3_code == eff_country & Time == 2000) %>% pull(PopTotal)
    
   pop <- pop * pre_2000_pop / pop_2000
   
  }

  area_per_1km2 <- as.numeric(sf::st_area(total_shp[total_shp$lctn_pr==lp_id,]$geometry)/1000/1000)
  
  rural_urban <- ifelse(
    pop/area_per_1km2 >= 1000,
    "urban",'rural'
  )
  
  pre_data$pop[row_idx] = pop
  pre_data$area_per_1km2[row_idx] = area_per_1km2
  pre_data$rural_urban[row_idx] = rural_urban
  
  cat(glue::glue("{year_tmp} is done for {lp_id}"))
}

processed_data <- read.csv("Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv") %>% 
  mutate(across(c("eff_tl", "eff_tr"), ~ char_to_date(.) %>% as.Date())) %>%
  mutate(study_year_tl = lubridate::year(eff_tl))%>% 
  left_join(
    pre_data, by = c('eff_country',"study_year_tl","lctn_pr")
  ) %>% 
  select(!study_year_tl)
#write.csv(processed_data,"Data/extracted_data/processed_cholera_climate_systematic_review_extract_into_sheets.csv",row.names = F)
