# cholera_climate_review
Systematic review of cholera and climate
## About
This project is a systematic review about assications between cholera and climate with a focus on precipitation, temperature, floods, and droughts.
## Overview
This repository includes raw data, generated datasets, analysis scripts, and summary notebooks.
## Data Description
### Climate zone raster data (`climatic_zone_raster`)
This folder contains a high-resolution global Koppen-Geiger climate classification maps available from https://www.gloh2o.org/koppen/ and a dictionary file to link the Koppen-Geiger climate classification to a more coarse classification. 
### Endemicity data (`endemicity`)
This folder contains a csv file with annual cholera cases from WHO Annual Reports and a csv file with the country ISO code and unique location identifier. Both data are avaiable from Cholera Taxonomy Databsae (https://cholera-taxonomy.middle-distance.com/)
### Extracted analytic datasets (`extracted_data`)
This folder contains the raw and processed datasets extracted from reviewing studies. 
### Human Development Index (HDI) data (`hdi_raster`)
This folder contains shapefile data from both the Cholera Taxonomy Database (https://cholera-taxonomy.middle-distance.com/) and The Humanitarian Data Exchange website (https://data.humdata.org/) and subnational HDI raster avaiable from UNDP GeoHub (https://geohub.data.undp.org/data/bb81caa483a0856029f980d2d7e89177)
### Inland and Coastal index data (`inland_coast`)
This folder contains the raw coastlines data from OpenStreetMap (https://osmdata.openstreetmap.de/data/coastlines.html) and processed distance to water data.
### Shapefile data (`shapefiles`)
This folder contains raw shapefile data for extracted locations and countires from both the Cholera Taxonomy Database (https://cholera-taxonomy.middle-distance.com/) and The Humanitarian Data Exchange website (https://data.humdata.org/), and the processed intersected shapefile data overlapping with the coastline data. 
