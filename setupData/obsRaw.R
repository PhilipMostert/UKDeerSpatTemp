library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(inlabru)
library(sf)
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)

species <- c("Capreolus capreolus", "Muntiacus reevesi", "Cervus elaphus", "Dama dama")

##Transects 1:5 T1
##Transects 6:10 T2
bto_obs <- read.csv('Data/Raw/4873 MOSTERT/4873_MOSTERT_BBS_TRANSECT_OBS.csv')

##Observations before 2014
bto_obs0 <- bto_obs %>%
  filter(YEAR > 2004, TRANSECT %in% "W")

#Transect 1
bto_obs1 <- bto_obs %>%
  filter(YEAR > 2004, TRANSECT %in% 1:5, DURATION_HRS < 2.5, DURATION_HRS > 0) %>%
  group_by(GRIDREF, YEAR, VISIT) %>%
  mutate(num_pres = n(),
         time_totalx = sum(DURATION_HRS/num_pres),
         n_transects = length(unique(TRANSECT)))

#Transect 2
bto_obs2 <- bto_obs %>%
  filter(YEAR > 2004, TRANSECT %in% 6:10, DURATION_HRS < 2.5, DURATION_HRS > 0) %>%
  group_by(GRIDREF, YEAR, VISIT) %>%
  mutate(num_pres = n(),
         time_totalx = sum(DURATION_HRS/num_pres),
         n_transects = length(unique(TRANSECT)))

bto_obs <- bind_rows(bto_obs0, bto_obs1, bto_obs2) %>%
  group_by(GRIDREF, YEAR) %>%
  mutate(time_total = sum(time_totalx/num_pres),
         n_visits = length(unique(VISIT)))

Times <- readxl::read_xlsx('Data/Raw//4873 MOSTERT/4873_MOSTERT_BBS_TRANSECT_TIMES.xlsx')

Times <- Times %>%
  filter(YEAR > 2004) %>%
  na.omit() %>%
  mutate(T1_START_TIME = hm(sub("(.{2})", "\\1 ", T1_START_TIME)),
         T1_END_TIME  = hm(sub("(.{2})", "\\1 ", T1_END_TIME)),
         T2_START_TIME = hm(sub("(.{2})", "\\1 ", T2_START_TIME)),
         T2_END_TIME = hm(sub("(.{2})", "\\1 ", T2_END_TIME)),
         T1_Duration = as.numeric(T1_END_TIME - T1_START_TIME, 'hours'),
         T2_Duration = as.numeric(T2_END_TIME - T2_START_TIME, 'hours'),
         T_Total = T1_Duration + T2_Duration) %>%
  filter(T1_Duration > 0, T2_Duration > 0,
         T1_Duration < 2.5, T2_Duration < 2.5) %>%
  group_by(GRIDREF, YEAR) %>%
  mutate(time_total = sum(T_Total),
         DURATION_HRS = 0,
         n_transects = 5,
         n_visits = length(unique(VISIT))) %>%
  select(GRIDREF, YEAR, VISIT, time_total, DURATION_HRS, n_visits, n_transects)


bto_comb <- rows_patch(x = bto_obs, y = Times, by = c('GRIDREF', "YEAR", 'VISIT'), unmatched = 'ignore') %>%
  as_tibble() %>%
  rename(species = ENGLISH_NAME) %>%
  filter(time_total > 0.08333333) %>%
  complete(GRIDREF, YEAR) %>%
  replace_na(list(time_total = 0, species = 'N')) %>%
  mutate(easting = rnrfa::osg_parse(GRIDREF)$easting,
         northing = rnrfa::osg_parse(GRIDREF)$northing,
         species = replace(species, species == "Roe Deer", "Capreolus capreolus"),
         species = replace(species, species == "Red Deer", "Cervus elaphus"),
         species = replace(species, species == "Fallow Deer", "Dama dama"),
         species = replace(species, species == "Reeves's Muntjac", "Muntiacus reevesi"))

bto_dat <- list()
bto_cv <- list()

for (spec in species) {

  bto_dat[[spec]] <- bto_comb %>%
    rename_all(tolower) %>%
    filter(!signs %in% c('L', 'S', 'L,S', 'F', 'F,L', 'F,S', 'F,L,S')) %>%
    filter(year >= 2005) %>%
    group_by(gridref, year, visit, easting, northing, time_total, n_visits) %>%
    summarize(count =  sum(`species` == spec), pres = sum(count > 0), species = spec) %>%
    ungroup() %>%
    group_by(year) %>%
    mutate(transects_year = sum(n_visits, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(gridref, year) %>%
    slice(1) %>%
    as.data.frame() %>%
    st_as_sf(coords = c('easting', 'northing'), crs = "epsg:27700") %>%
    st_transform(crs) %>%
    st_filter(UKsimp) %>%
    select(gridref, year, pres, transects_year, species, n_visits)

  bto_cv[[spec]] <- bto_comb %>%
    rename_all(tolower) %>%
    filter(signs %in% c('F', 'F,L', 'F,S', 'F,L,S')) %>%
    filter(year >= 2005 & year <= 2024) %>%
    group_by(gridref, year, visit, easting, northing, time_total, n_visits) %>%
    summarize(count =  sum(`species` == spec), pres = sum(count > 0), species = spec) %>%
    ungroup() %>%
    group_by(gridref, year) %>%
    slice(1) %>%
    filter(count > 0) %>%
    as.data.frame() %>%
    st_as_sf(coords = c('easting', 'northing'), crs = "epsg:27700") %>%
    st_transform(crs) %>%
    st_filter(UKsimp) %>%
    select(gridref, year, pres, species)

}

saveRDS(do.call(rbind, bto_dat), 'Data/BTO.rds')
saveRDS(do.call(rbind, bto_cv), 'Data/BTOCV.rds')

iNat <- list()

for (spec in species) {

  rawObs <- rgbif::occ_data(datasetKey = '50c9509d-22c7-4a22-a47d-8c48425ef4a7',
                         scientificName = spec, limit = 20000,
                         coordinateUncertaintyInMeters = '0,1000', #'0,1000'
                         country = 'GB')

  iNat[[spec]] <- rawObs$data

}

iNat <- lapply(iNat, function(x){

  y <- st_as_sf(x, coords = c('decimalLongitude','decimalLatitude'), crs = 'WGS84') %>%
    st_transform(crs) %>%
    st_filter(UKsimp) %>%
    filter(year > 1994 & year <= 2025) %>%
    select(year, species)
  y

})

saveRDS(do.call(rbind, iNat), 'Data/iNat.rds')


Mam <- read.csv('Data/Raw/MAMATLAS//MAMATLAS.csv')
MamObs <- list()
for (spec in species) {

  MamObs[[spec]] <- Mam %>%
    filter(`State.Province` %in% c("England", "Scotland", "Wales")) %>%
    rename(year = `Start.date.year`, species = `Scientific.name`) %>%
    filter(species == spec) %>%
    filter(year > 1994 & year <= 2025) %>%
    filter(Coordinate.uncertainty..m. <= 1000) %>% #New
    st_as_sf(coords = c('Longitude..WGS84.', 'Latitude..WGS84.'), crs = 'WGS84') %>%
    st_transform(crs) %>%
    st_filter(UKsimp) %>%
    select(year, species)

}
saveRDS(do.call(rbind, MamObs), 'Data/Mammal.rds')
