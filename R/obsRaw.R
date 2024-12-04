
library(dplyr)
library(tidyr)
library(ggplot2)
library(inlabru)
library(sf)
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)

species <- c("Capreolus capreolus", "Muntiacus reevesi", "Cervus elaphus", "Dama dama")

btoObs <- list()
bto_data <- read.csv('Data/Raw/AvianNew/records-2024-10-07.csv')
for (spec in species) {

  btoObs[[spec]] <- bto_data %>%
    filter(`State.Province` %in% c("England", "Scotland", "Wales")) %>%
    rename(OSGR10K = `OSGR.10km`, year = `Start.date.year`, species = `Scientific.name`) %>%
    group_by(OSGR10K, year) %>%
    summarise(listLength = length(unique(`species`)),
              numPres = sum(`species` == spec),
              pres = ifelse(numPres >0, 1, 0),
              obsTotal = sum(`Start.date.month` %in% 1:12),
              species = species, year = year) %>% slice(1) %>% ungroup() %>%#.groups = "drop"
    mutate(#presence = ifelse(presence > 0, 1, 0),
      species =  spec,
      lon = rnrfa::osg_parse(OSGR10K, coord_system = 'WGS84')$lon,
      lat = rnrfa::osg_parse(OSGR10K, coord_system = 'WGS84')$lat) %>%
    st_as_sf(coords = c('lon', 'lat'), crs = 'WGS84') %>%
    st_transform(crs) %>%
    st_filter(UKsimp)

}
saveRDS(do.call(rbind, btoObs), 'Data/BTO.rds')

iNat <- list()

#iNat 50c9509d-22c7-4a22-a47d-8c48425ef4a7
#Mam 67e2335c-4303-4ee3-9bfe-ecfe65662629
#MamWeb 8eafd1af-58c2-49cf-b093-6796b9b38c8e

for (spec in species) {

  rawObs <- rgbif::occ_data(datasetKey = '50c9509d-22c7-4a22-a47d-8c48425ef4a7',
                         scientificName = spec, limit = 20000, country = 'GB',
                         coordinateUncertaintyInMeters = '0,1000')

  iNat[[spec]] <- rawObs$data

}

iNat <- lapply(iNat, function(x){

  y <- st_as_sf(x, coords = c('decimalLongitude','decimalLatitude'), crs = 'WGS84') %>%
    st_transform(crs) %>%
    st_filter(UKsimp) %>%
    filter(year > 1994 & year < 2022) %>%
    select(year, species)
  y

})

saveRDS(do.call(rbind, iNat), 'Data/iNat.rds')


Mam <- read.csv('Data/Raw/MammalAtlasNew/MammalAtlas.csv')
MamObs <- list()
for (spec in species) {

  MamObs[[spec]] <- Mam %>%
    filter(`State.Province` %in% c("England", "Scotland", "Wales")) %>%
    rename(year = `Start.date.year`, species = `Scientific.name`) %>%
    filter(species == spec) %>%
    filter(year > 1994 & year < 2022) %>%
    st_as_sf(coords = c('Longitude..WGS84.', 'Latitude..WGS84.'), crs = 'WGS84') %>%
    st_transform(crs) %>%
    st_filter(UKsimp) %>%
    select(year, species)

}
saveRDS(do.call(rbind, MamObs), 'Data/Mammal.rds')

MamWeb <- read.csv('Data/Raw/MammalWebNew/MammalWeb.csv')
MamWebObs <- list()
for (spec in species) {

  MamWebObs[[spec]] <- MamWeb %>%
    filter(`State.Province` %in% c("England", "Scotland", "Wales")) %>%
    rename(year = `Start.date.year`, species = `Scientific.name`) %>%
    filter(species == spec) %>%
    filter(year > 1994 & year < 2022) %>%
    st_as_sf(coords = c('Longitude..WGS84.', 'Latitude..WGS84.'), crs = 'WGS84') %>%
    st_transform(crs) %>%
    st_filter(UKsimp) %>%
    select(year, species)

}
saveRDS(do.call(rbind, MamWebObs), 'Data/MammalWeb.rds')

#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2664.13622
GBC <- read.csv('Data/Raw/Noble_etal_2012_DeerDensityMaps2005-09.csv')


extract_mean <- function(x) {

  if (x == 'Not reported') numbers <- NA#0
  else numbers <- as.numeric(unlist(strsplit(gsub("[^0-9]+", " ", x), " ")))
  mean(numbers)
}

GBC_final <- GBC %>%
  rowwise() %>%
  mutate(across(2:5, extract_mean),
         easting = rnrfa::osg_parse(os_grid)$easting,
         northing = rnrfa::osg_parse(os_grid)$northing) %>%
  ungroup() %>%
  pivot_longer(cols = 2:5, names_to = "species", values_to = "NGC") %>%
  drop_na() %>%
  mutate(species = replace(species, species == "roe_deer", "Capreolus capreolus"),
         species = replace(species, species == "red_deer", "Cervus elaphus"),
         species = replace(species, species == "fallow_deer", "Dama dama"),
         species = replace(species, species == "muntjac", "Muntiacus reevesi"),
         lon = rnrfa::osg_parse(os_grid, 'WGS84')$lon,
         lat = rnrfa::osg_parse(os_grid, 'WGS84')$lat) %>%
  sf:::st_as_sf(coords = c('lon', 'lat'), crs  ='WGS84') %>%
  st_transform(crs = crs) %>%
  st_filter(UKsimp)
saveRDS(GBC_final, 'Data/GBC.rds')

uniqueCoords <- GBC_final %>%
  st_coordinates() %>%
  unique() %>%
  data.frame() %>%
  st_as_sf(coords = c('X', 'Y'), crs = crs)

GBCPred <- list()
for (spec in species) {

  uniqueCoords$NGC <- 0

  uniqueCoords[unlist(st_equals(GBC_final[GBC_final$species == spec,], y = uniqueCoords)),]$NGC <- GBC_final[GBC_final$species == spec,]$NGC
  uniqueCoords$species <- spec
  GBCPred[[spec]] <- uniqueCoords

}
saveRDS(do.call(rbind, GBCPred), 'Data/GBCPredict.rds')


#dat <- read.csv(file = 'Data/Raw/AvianNew/Avian3.csv')

#dt <- st_as_sf(dat, coords = c('decimalLongitude.processed','decimalLatitude.processed'), crs = 'WGS84')

#dt <- dt %>% rename(species = scientificName.processed,
#                    year = year.processed) %>%
#  filter(year > 1994 & year < 2022)

#unique_species_locations <- dt %>%
#  select(species, geometry) %>%
#  distinct()

#all_time_periods <- unique(dt$year)

#species_time_grid <- expand.grid(species = unique(unique_species_locations$species),
#                                 year = all_time_periods,
#                                 stringsAsFactors = FALSE)

#complete_grid <- species_time_grid %>%
#  left_join(unique_species_locations, by = "species")# %>%
  #mutate(presence = NA)

##presence_data <- dt %>%
#  mutate(presence = 1) %>%
#  select(species, year, presence, 'geometry')

#complete_data <- complete_grid %>%
#  left_join(presence_data, by = c("species", "year", 'geometry')) %>%
#  mutate(presence = ifelse(is.na(presence), 0, presence))  %>%
#  st_as_sf()

#saveRDS(complete_data, 'Data/Mammal.rds')
