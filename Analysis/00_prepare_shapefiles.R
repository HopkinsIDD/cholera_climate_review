# This script is used to generate shapefiles linked to extracted data ----

# Load R packages ----
library(tidyverse)
library(stringdist)

# Import extracted spreadsheet ----
extracted_data <- readxl::read_excel('~cholera_climate_systematic_review_extract_into_sheets.xlsx', sheet = "estimates_0329") %>% 
  mutate(
    row_uid = 1:n()
  )

# Remove locations with regions/states/provinces and multi-locations since they are unable to be linked to shapefiles ----
### 71 rows from 9 studies
composite_loc <- extracted_data[grepl("(region|province|state|district|commune|multi|country)", extracted_data$eff_loc, ignore.case = TRUE), ]
#write.csv(composite_loc,"filtered_unknown_composite_loc.csv",row.names = F)

# Filter out rows with single locations ----
filtered_data <- extracted_data[!grepl("(region|province|state|district|commune|multi|country)", extracted_data$eff_loc, ignore.case = TRUE), ] %>% 
  rowwise() %>% 
  mutate(std_loc = paste(eff_country, eff_loc, sep = "::"),
         std_loc = gsub("::NA", "", std_loc))  # Remove any "NA" text

# Load location data from the outbreak study ----
loc_tax <- read.csv("location.csv")

# Create a function to find the best match for each location in filtered_data ----
find_best_match <- function(eff_loc, loc_tax_names) {
  # Compute similarity scores for each qualified_name
  similarity_scores <- stringdist::stringdist(eff_loc, loc_tax_names, method = "jw")
  # Find the index of the smallest distance
  best_match_index <- which.min(similarity_scores)
  # Return the best match
  return(loc_tax_names[best_match_index])
}
## Apply the function to match each eff_loc in filtered_data to loc_tax ----
filtered_data$best_match <- sapply(filtered_data$std_loc, find_best_match, loc_tax_names = loc_tax$qualified_name)
#write.csv(filtered_data,"filtered_best_match.csv",row.names=F)
filtered_data <- read.csv("filtered_best_match.csv")

# Load location period data from the outbreak study ----
lp <- read.csv("location_periods.csv") %>% 
  mutate(id = as.character(id))
filtered_data$location_id = NA
filtered_data$location_period_id = NA

for (idx in 1:nrow(filtered_data)) {
  filtered_data$location_id[idx] = loc_tax[loc_tax$qualified_name == filtered_data$best_match[idx],]$id
  
  if(nrow(lp[lp$location_id == filtered_data$location_id[idx],])!=0){
    filtered_data$location_period_id[idx] = lp[lp$location_id == filtered_data$location_id[idx],]$id
  }
}

# Extract shapefiles from cholera taxonomy database ----
## non-composite locations (specify the status by manually checking whether the best_match looks valid) ----
### done means that it looks good
### missing_loc means that it needs location audit
### missing_shp means that it needs shapefiles
data_with_shps <- read.csv("filtered_best_match.csv") %>% 
  filter(status %in% c("done"))
lp_row <- data_with_shps %>% 
  distinct(row_uid,location_period_id)

## composite locations (there are still some composite locations in the filtered data) ----
data_composite_loc <- read.csv("filtered_best_match.csv") %>% 
  filter(status %in% c("composite_loc"))

### Create a unique numeric identifier for each composite location ----
data_composite_loc$composite_location_period_id <- paste("composite_loc",as.numeric(factor(data_composite_loc$location_period_id)),sep="-")

### Create a list of child location identifiers for each composite location ----
data_composite_loc$child_locations <- sapply(data_composite_loc$location_period_id, function(x) {
  if (grepl("\\|", x)) {
    child_ids <- unlist(strsplit(x, "\\|"))
    child_ids
  } else {
    return(NA)  # No composite location, no child
  }
})


lp_row_composite <- data_composite_loc %>% 
  distinct(row_uid,composite_location_period_id,child_locations)

lp_row_composite_split <- lp_row_composite %>%
  unnest(child_locations) %>%
  arrange(row_uid) %>% # Arrange by row_uid (optional)
  mutate(location_period_id = child_locations)

### location period ids that need to be pulled from the database ----
lp_row <- rbind(lp_row, lp_row_composite_split%>% distinct(row_uid,location_period_id))

#write.csv(lp_row,"lp_from_database.csv",row.names = F)

## Get shapefiles on idm ----
lp_row <- read.csv("lp_from_database.csv")
climate_shp <- data.frame()
for (tmp_lp in unique(lp_row$location_period_id[!lp_row$location_period_id == 18105])){
  tmp_shp <- get_shp(
    location_period_id = tmp_lp,
    username = username,
    password = password,
    host="db.cholera-taxonomy.middle-distance.com",
    port=5432,
    dbname='CholeraTaxonomy_production'
  ) %>% 
    select(lctn_pr,geometry)
  
  if(length(climate_shp)==0){
    climate_shp=tmp_shp
  } else {
    climate_shp=rbind(climate_shp,tmp_shp)    
  }
  print(tmp_lp)
}
#saveRDS(climate_shp, "climate_shp.rds")

## Combine the shapefiles for composite locations ----
lp_row <- read.csv("lp_from_database.csv") 
climate_shp <- readRDS("climate_shp.rds")

# Join the datasets based on location_period_id == lctn_pr
merged_data <- lp_row %>%
  left_join(climate_shp, by = c("location_period_id" = "lctn_pr"))

# Union geometries for each row_uid
unioned_geometries <- merged_data %>%
  group_by(row_uid) %>%
  summarise(geometry = st_union(geometry), .groups = 'drop')
#st_write(unioned_geometries, "unioned_climate_shapefile.shp")

filtered_data$new_location_period_id = NA
final_climate_shp <- data.frame()
for (row_idx in unique(unioned_geometries$row_uid)) {
  if(row_idx %in% unique(lp_row_composite$row_uid)){
    lctn_pr = lp_row_composite$composite_location_period_id[lp_row_composite$row_uid==row_idx]
    geometry = unioned_geometries$geometry[unioned_geometries$row_uid==row_idx]
  } else{
    lctn_pr = lp_row$location_period_id[lp_row$row_uid==row_idx]
    geometry = unioned_geometries$geometry[unioned_geometries$row_uid==row_idx]
  }
  
  filtered_data$new_location_period_id[filtered_data$row_uid==row_idx] = lctn_pr
  tmp_shp = data.frame(lctn_pr = lctn_pr,geometry=geometry) %>% sf::st_as_sf()
  
  if(nrow(final_climate_shp)==0){
    final_climate_shp=tmp_shp
  }else{
    final_climate_shp = rbind(final_climate_shp,tmp_shp)
  }
  print(row_idx)
}

## Link location period ids to the extracted data ----
filtered_data$new_location_period_id = NA
for (row_idx in unique(unioned_geometries$row_uid)) {
  if(row_idx %in% unique(lp_row_composite$row_uid)){
    lctn_pr = lp_row_composite$composite_location_period_id[lp_row_composite$row_uid==row_idx]
  } else{
    lctn_pr = lp_row$location_period_id[lp_row$row_uid==row_idx]
  }  
  filtered_data$new_location_period_id[filtered_data$row_uid==row_idx] = lctn_pr
  print(row_idx)
}

final_climate_shp1 <- final_climate_shp %>% 
  group_by(lctn_pr) %>% 
  slice(1)
#saveRDS(final_climate_shp1,"climate_unique_shp.rds")

#write.csv(filtered_data,"filtered_data_shp.csv",row.names=F)

total_shp <- readRDS("climate_unique_shp.rds")

# Look for shapefiles from other resources for those locations without shapefiles ----
## bgd matlab
bgd_matlab <- sf::st_read("new_shps/bgd_admbnda_adm3_bbs_20201113.shp") %>% 
  filter(ADM3_EN %in% c("Matlab Dakshin","Matlab Uttar")) %>% 
  summarise(
    geometry = st_union(geometry)
  ) %>% mutate(
    lctn_pr = "new-1"
  )

final_climate_shp2 <- rbind(final_climate_shp1,bgd_matlab)

## malaysia
mys_Terengganu <- sf::st_read("new_shps/mys_admbnda_adm1_unhcr_20210211.shp") %>% 
  filter(ADM1_EN == "Terengganu") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-2"
  )

mys_Sarawak <- sf::st_read("new_shps/mys_admbnda_adm1_unhcr_20210211.shp") %>% 
  filter(ADM1_EN == "Sarawak") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-3"
  )

mys_Kedah <- sf::st_read("new_shps/mys_admbnda_adm1_unhcr_20210211.shp") %>% 
  filter(ADM1_EN == "Kedah") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-4"
  )

mys_Kelantan <- sf::st_read("new_shps/mys_admbnda_adm1_unhcr_20210211.shp") %>% 
  filter(ADM1_EN == "Kelantan") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-5"
  )

mys_Selangor <- sf::st_read("new_shps/mys_admbnda_adm1_unhcr_20210211.shp") %>% 
  filter(ADM1_EN == "Selangor") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-6"
  )

final_climate_shp3 <- rbind(final_climate_shp2,mys_Terengganu,mys_Sarawak,mys_Kedah,mys_Kelantan,mys_Selangor)
#saveRDS(final_climate_shp3,"climate_unique_shp3.rds")

## peru
per_Castilla <- sf::st_read("new_shps/per_admbnda_adm2_ign_20200714.shp") %>% 
  filter(ADM2_ES == "Castilla") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-7"
  )


per_Catacaos <- sf::st_read("new_shps/per_admbnda_adm3_ign_20200714.shp") %>% 
  filter(ADM3_ES == "Catacaos") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-8"
  )

per_Chulucanas <- sf::st_read("new_shps/per_admbnda_adm3_ign_20200714.shp") %>% 
  filter(ADM3_ES == "Chulucanas") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-9"
  )



per_LaArena <- sf::st_read("new_shps/per_admbnda_adm3_ign_20200714.shp") %>% 
  filter(ADM3_ES == "La Arena") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-10"
  )


per_LaUnion <- sf::st_read("new_shps/per_admbnda_adm2_ign_20200714.shp") %>% 
  filter(ADM2_ES == "La Union") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-11"
  )

per_Paita <- sf::st_read("new_shps/per_admbnda_adm2_ign_20200714.shp") %>% 
  filter(ADM2_ES == "Paita") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-12"
  )


per_Piura <- sf::st_read("new_shps/per_admbnda_adm2_ign_20200714.shp") %>% 
  filter(ADM2_ES == "Piura") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-13"
  )


per_RinconadaLlicuar <- sf::st_read("new_shps/per_admbnda_adm3_ign_20200714.shp") %>% 
  filter(ADM3_ES == "Rinconada Llicuar") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-14"
  )


per_Sechura <- sf::st_read("new_shps/per_admbnda_adm3_ign_20200714.shp") %>% 
  filter(ADM3_ES == "Sechura") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-15"
  )

final_climate_shp4 <- rbind(final_climate_shp3,per_Castilla,per_Catacaos,per_Chulucanas,per_LaArena,per_LaUnion,per_Paita,per_Piura,per_RinconadaLlicuar,per_Sechura)
#saveRDS(final_climate_shp4,"climate_unique_shp4.rds")

## iran
irn_Chabahar <- sf::st_read("new_shps/irn_admbnda_adm2_unhcr_20190514.shp") %>% 
  filter(ADM2_EN == "Chabahar") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-16"
  )
final_climate_shp5 <- rbind(final_climate_shp4,irn_Chabahar)
#saveRDS(final_climate_shp5,"climate_unique_shp5.rds")

## vitnam
vet_Hue <- sf::st_read("new_shps/vnm_admbnda_adm2_gov_20201027.shp") %>% 
  filter(ADM2_EN == "Hue") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-17"
  )

vet_NhaTrang <- sf::st_read("new_shps/vnm_admbnda_adm2_gov_20201027.shp") %>% 
  filter(ADM2_EN == "Nha Trang") %>% 
  select(geometry) %>%
  mutate(
    lctn_pr = "new-18"
  )

final_climate_shp6 <- rbind(final_climate_shp5,vet_Hue,vet_NhaTrang)
#saveRDS(final_climate_shp6,"climate_unique_shp6.rds")

## nigeria
nga_ibadan <- sf::st_read("new_shps/new_lga_nigeria_2003.shp") %>% 
  filter(LGA %in% c("Ibadan North","Ibadan North-East","Ibadan North-West","Ibadan South-East","Ibadan South-West")) %>% 
  summarise(geometry = sf::st_union(geometry)) %>% 
  mutate(
    lctn_pr = "new-19"
  )

final_climate_shp7 <- rbind(final_climate_shp6,nga_ibadan)
#saveRDS(final_climate_shp7,"climate_unique_shp7.rds")

## india
ind_kolkata <- sf::st_read("new_shps/geoBoundaries-IND-ADM2.geojson") %>% 
  filter(shapeName == "Kolkata") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-20"
  )
final_climate_shp8 <- rbind(final_climate_shp7,ind_kolkata)
#saveRDS(final_climate_shp8,"climate_unique_shp8.rds")

## China
chn_zhejiang <- sf::st_read("new_shps/chn_admbnda_adm1_ocha_2020.shp") %>% 
  filter(ADM1_EN == "Zhejiang Province") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-21"
  )

final_climate_shp9 <- rbind(final_climate_shp8,chn_zhejiang)
#saveRDS(final_climate_shp9,"climate_unique_shp9.rds")

## India Vellore
ind_Vellore <- sf::st_read("new_shps/geoBoundaries-IND-ADM2.geojson") %>% 
  filter(shapeName == "Vellore")%>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-22"
  )
final_climate_shp10 <- rbind(final_climate_shp9,ind_Vellore)
#saveRDS(final_climate_shp10,"climate_unique_shp10.rds")

ind_Ludhiana <- sf::st_read("new_shps/geoBoundaries-IND-ADM2.geojson") %>% 
  filter(shapeName == "Ludhiana")%>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-23"
  )

ind_Sirsa <- sf::st_read("new_shps/geoBoundaries-IND-ADM2.geojson") %>% 
  filter(shapeName == "Sirsa")%>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-24"
  )

final_climate_shp11 <- rbind(final_climate_shp10,ind_Ludhiana,ind_Sirsa)
#saveRDS(final_climate_shp11,"climate_unique_shp11.rds")

## Leckebusch2015: Kano, Sokoto and Gusau
nga_united <- sf::st_read("new_shps/new_lga_nigeria_2003.shp") %>% 
  filter(STATE %in% c("Kano","Sokoto") | LGA =="Gusau") %>% 
  summarise(geometry = sf::st_union(geometry)) %>% 
  mutate(
    lctn_pr = "new-25"
  )

plot(nga_united)
final_climate_shp12 <- rbind(final_climate_shp11,nga_united)
#saveRDS(final_climate_shp12,"climate_unique_shp12.rds")

## Jutla2013: ind and pak: across 9 regions
#(delhi, lahore, ludhiana, sialkot, rawalpindi, peshawar, dera Ismail Khan, multan, and sirsa)
#21252, 7446, "new-23", 6028,8353,4417,3563,3282,"new-24"

ind_pak <- final_climate_shp12 %>% 
  filter(
    lctn_pr %in% c(21252, 7446, "new-23", 6028,8353,4417,3563,3282,"new-24")
  ) 


ind_pak_shp<-sf::st_union(ind_pak) %>% 
  sf::st_as_sf()%>%
  rename(geometry=x) %>%
  mutate(lctn_pr = "new-26")

final_climate_shp13 <- rbind(final_climate_shp12,ind_pak_shp)
#saveRDS(final_climate_shp13,"climate_unique_shp13.rds")

## Oloukoi2014:
nga_united2 <- sf::st_read("new_shps/new_lga_nigeria_2003.shp") %>% 
  filter(LGA %in% c("Iseyin", "Kajola","Saki West")) %>% 
  summarise(geometry = sf::st_union(geometry)) %>% 
  mutate(
    lctn_pr = "new-27"
  )

final_climate_shp14 <- rbind(final_climate_shp13,nga_united2)
#saveRDS(final_climate_shp14,"climate_unique_shp14.rds")

## Yang2022
### Difference climatic zones in China
chn_temperate<- sf::st_read("new_shps/chn_admbnda_adm1_ocha_2020.shp") %>% 
  filter(ADM1_EN %in% c(
    "Xinjiang Uygur Autonomous Region", 
    "Gansu province",
    "Inner Mongolia Autonomous Region",
    "Heilongjiang Province",
    "Jilin Province",
    "Liaoning Province",
    "Beijing Municipality",
    "Tianjin Municipality",
    "Shandong Province",
    "Hebei Province",
    "Shanxi Province",
    "Henan Province",
    "Shaanxi Province",
    "Ningxia Hui Autonomous Region")) %>% 
  select(geometry)

chn_temperate_union <- sf::st_union(chn_temperate$geometry) %>% 
  sf::st_as_sf() %>% 
  rename(geometry =x )%>%
  mutate(
    lctn_pr = "new-28"
  )


chn_subtropical<- sf::st_read("new_shps/chn_admbnda_adm1_ocha_2020.shp") %>% 
  filter(ADM1_EN %in% c(
    "Sichuan Province",
    "Yunnan Province",
    "Guizhou Province",
    "Chongqing Municipality",
    "Hubei Province",
    "Hunan Province",
    "Guangxi Zhuang Autonomous Region",
    "Anhui Province",
    "Jiangxi Province",
    "Guangdong Province",
    "Jiangsu Province",
    "Shanghai Municipality",
    "Zhejiang Province",
    "Fujian Province",
    "Guangdong Province",
    "Hong Kong Special Administrative Region",
    "Macao Special Administrative Region")) %>% 
  select(geometry) 

chn_subtropical_union <- sf::st_union(chn_subtropical$geometry) %>% 
  sf::st_as_sf() %>% 
  rename(geometry =x )%>%
  mutate(
    lctn_pr = "new-29"
  )

chn_tropical<- sf::st_read("new_shps/chn_admbnda_adm1_ocha_2020.shp") %>% 
  filter(ADM1_EN %in% c(
    "Hainan Province")) %>% 
  select(geometry) 

chn_tropical_union <- chn_tropical %>%
  mutate(
    lctn_pr = "new-30"
  )

chn_alpine<- sf::st_read("new_shps/chn_admbnda_adm1_ocha_2020.shp") %>% 
  filter(ADM1_EN %in% c(
    "Tibet Autonomous Region","Qinghai Province")) %>% 
  select(geometry) 

chn_alpine_union <-sf::st_union(chn_alpine$geometry) %>% 
  sf::st_as_sf() %>% 
  rename(geometry =x )%>%
  mutate(
    lctn_pr = "new-31"
  )

final_climate_shp15 <- rbind(final_climate_shp14,chn_temperate_union,chn_subtropical_union,chn_tropical_union,chn_alpine_union)
#saveRDS(final_climate_shp15,"climate_unique_shp15.rds")

## country -level shpaefiles:
# uga, gha, yem, tza, nga, china
uga<-rgeoboundaries::gb_adm0(country = "UGA") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-32"
  )

gha<-rgeoboundaries::gb_adm0(country = "GHA") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-33"
  )

yem<-rgeoboundaries::gb_adm0(country = "YEM") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-34"
  )

tza<-rgeoboundaries::gb_adm0(country = "TZA") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-35"
  )

nga<-rgeoboundaries::gb_adm0(country = "NGA") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-36"
  )


chn<-rgeoboundaries::gb_adm0(country = "CHN") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-37"
  )

final_climate_shp16 <- rbind(final_climate_shp15,uga,gha,yem,tza,nga,chn)
#saveRDS(final_climate_shp16,"climate_unique_shp16.rds")

hti<-rgeoboundaries::gb_adm0(country = "HTI") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-38"
  )

final_climate_shp17 <- rbind(final_climate_shp16,hti)
#saveRDS(final_climate_shp17,"climate_unique_shp17.rds")

## vietnam
vnm<-rgeoboundaries::gb_adm0(country = "VNM") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-39"
  )

final_climate_shp18 <- rbind(final_climate_shp17,vnm)
#saveRDS(final_climate_shp18,"climate_unique_shp18.rds")

## kenya
ken<-rgeoboundaries::gb_adm0(country = "KEN") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-40"
  )

final_climate_shp19 <- rbind(final_climate_shp18,ken)
#saveRDS(final_climate_shp19,"climate_unique_shp19.rds")

## cameroon different climatic zones
cmr <- sf::st_read("new_shps/cmr_admbnda_adm1_inc_20180104.shp")

### Equitorial Monsoon: Northwest, West, Southwest, and Littoral
cmr_em <- sf::st_read("new_shps/cmr_admbnda_adm1_inc_20180104.shp") %>% 
  filter(ADM1_EN %in% c("North-West","West","South-West","Littoral")) %>% 
  select(geometry) 

cmr_em_shp <- sf::st_union(cmr_em$geometry) %>% 
  sf::st_as_sf() %>%
  rename(geometry = x) %>% 
  mutate(
    lctn_pr = "new-41"
  )

### Guinea Equatorial: Centre, East, and South
cmr_ge <- sf::st_read("new_shps/cmr_admbnda_adm1_inc_20180104.shp") %>% 
  filter(ADM1_EN %in% c("Centre","East","South")) %>% 
  select(geometry) 

cmr_ge_shp <- sf::st_union(cmr_ge$geometry) %>% 
  sf::st_as_sf() %>%
  rename(geometry = x) %>% 
  mutate(
    lctn_pr = "new-42"
  )

### Sudano-Sahelian: far-north
cmr_fn <- sf::st_read("new_shps/cmr_admbnda_adm1_inc_20180104.shp") %>% 
  filter(ADM1_EN %in% c("Far-North")) %>% 
  select(geometry) 

cmr_fn_shp <- sf::st_union(cmr_fn$geometry) %>% 
  sf::st_as_sf() %>%
  rename(geometry = x) %>% 
  mutate(
    lctn_pr = "new-43"
  )

### Tropical Humid: North and Adamawa
cmr_th <- sf::st_read("new_shps/cmr_admbnda_adm1_inc_20180104.shp") %>% 
  filter(ADM1_EN %in% c("North","Adamawa")) %>% 
  select(geometry) 

cmr_th_shp <- sf::st_union(cmr_th$geometry) %>% 
  sf::st_as_sf() %>%
  rename(geometry = x) %>% 
  mutate(
    lctn_pr = "new-44"
  )

final_climate_shp20 <- rbind(final_climate_shp19,cmr_em_shp,cmr_ge_shp,cmr_fn_shp,cmr_th_shp)
#saveRDS(final_climate_shp20,"climate_unique_shp20.rds")

# African continent
afr <- rgeoboundaries::gb_adm0(country = c("bdi","dji","eri","eth","ken","mwi","moz","rwa","som","tza","uga","zmb","ago","cmr","caf","tcd","cog","cod","gnq","gab","bwa","LSO","nam","zaf","swz","zwe","ben","bfa","civ","gmb","gha","gin","gnb","lbr","mli","mrt","ner","nga","sen","sle","tgo"))
afr_shp <- sf::st_union(afr$geometry) %>% 
  sf::st_as_sf() %>%
  rename(geometry = x) %>% 
  mutate(
    lctn_pr = "new-45"
  )

final_climate_shp21 <- rbind(final_climate_shp20,afr_shp)
#saveRDS(final_climate_shp21,"climate_unique_shp21.rds")

# Rename climate_unique_shp21 to location_shapefiles.rds in Data/shapefiles folder in cholera_climate_review Github repo. ----
# HTI: Artibonite
hti_Artibonite <- rgeoboundaries::gb_adm1("HTI") %>% 
  filter(shapeName == "Artibonite Department") %>% 
  select(geometry) %>% 
  mutate(
    lctn_pr = "new-47"
  )
total_shps <- rbind(total_shps,hti_Artibonite)
#saveRDS(total_shps,"cholera_climate_review/Data/shapefiles/location_shapefiles.rds")
