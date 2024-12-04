library(sf)
library(terra)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(ggplot2)
library(inlabru)
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)

for (i in 1:26) {

file <- list.files(path = 'Climate/temp/')[i]
xx <- terra::rast(paste0('Climate/temp/', file))
xx <- project(xx, crs(UK))
xx <- mask(crop(xx, UK), UK)

if (i == 1) temp <- xx
else temp <- c(temp, xx)

}

names(temp) <- 1995:2020
temp <- project(temp, crs)

for (i in 1:26) {

  file <- list.files(path = 'Climate/prec/')[i]
  xx <- terra::rast(paste0('Climate/prec/', file))
  xx <- project(xx, crs(UK))
  xx <- mask(crop(xx, UK), UK)

  if (i == 1) prec <- xx
  else prec <- c(prec, xx)

}

names(prec) <- 1995:2020
prec <- project(prec, crs)

  lc2000 <- crop(terra::rast('Climate/landcover/2000/U2006_CLC2000_V2020_20u1.tif'), UK)
  lc2006 <- crop(terra::rast('Climate/landcover/2006/U2012_CLC2006_V2020_20u1.tif'), UK)
  lc2012 <- crop(terra::rast('Climate/landcover/2012/U2018_CLC2012_V2020_20u1_raster100m.tif'), UK)
  lc2018 <- crop(terra::rast('Climate/landcover/2018/U2018_CLC2018_V2020_20u1.tif'), UK)

  newValues <- c(rep('Urban', 11),
                 rep('Arable land', 6),
                 'Arable land', #APastures
                 'Arable land',
                 'Arable land',
                 'Arable land',
                 'Arable land',
                 'Forest', #Forest Broad-leaved forest
                 'Forest', #Forest Coniferous forest
                 NA, #Mixed forest
                 'Natural grasslands', #Natural grasslands
                 'Moors and heathland', # Moors and heathland
                 'Natural grasslands', # Sclerophyllous vegetation
                 'Natural grasslands', # Transitional woodland-shrub
                 rep(NA, 2), #Rocks_sands or NA
                 'Natural grasslands', #Sparsely vegetated areas
                 rep(NA, 2),
                 'Wetlands', #Wetlands #Inland marshes
                 'Wetlands',# Peat bogs
                 rep('Wetlands', 8),
                 NA) ##Wetlands

  values(lc2000) <- newValues[c(values(lc2000))]
  values(lc2006) <- newValues[c(values(lc2006))]
  values(lc2012) <- newValues[c(values(lc2012))]
  values(lc2018) <- newValues[c(values(lc2018))]
  lc <- c(lc2000, lc2006, lc2012, lc2018)
  lc <- aggregate(lc, 10, fun = 'modal', na.rm = TRUE)
  names(lc) <- c('lc2000', 'lc2006', 'lc2012', 'lc2018')
  lc <- mask(crop(lc, UK), UK)
  lc <- project(lc, crs)

elevation <- elevatr::get_elev_raster(locations = UK, z = 7)
elevation <- mask(crop(elevation, UK), UK)
elevation <- project(rast(elevation), crs)
elevation <- clamp(elevation, 0)
units(elevation) <- 'm'
names(elevation) <- 'elevation'

LC <- terra::rast('Climate/LC/data/LCM2015_GB_1km_percent_cover_aggregate_class.tif')

Arable <- project(LC$Layer_3, crs(UK))
Built <- project(LC$Layer_10, crs(UK))
Coniferous <- project(LC$Layer_2, crs(UK))
Grassland <- project(LC$Layer_4, crs(UK))

Arable <- mask(crop(Arable, UK), UK)
Arable <- project(Arable, crs)
names(Arable) <- 'Arable'

Coniferous <- mask(crop(Coniferous, UK), UK)
Coniferous <- project(Coniferous, crs)
names(Coniferous) <- 'Coniferous'

Grassland <- mask(crop(Grassland, UK), UK)
Grassland <- project(Grassland, crs)
names(Grassland) <- 'Grassland'

Built <- mask(crop(Built, UK), UK)
Built <- project(Built, crs)
names(Built) <- 'Built'


covs <- list(Temperature = temp,
             Precipitation = prec,
             LandCover = lc,
             Elevation = elevation)

#saveRDS(covs, 'Data/SpatRasters.rds')
covs2 <- readRDS('Data/SpatRasters.rds')
terra::writeRaster(temp, 'Data/temp.tiff', overwrite = TRUE)
terra::writeRaster(prec, 'Data/prec.tiff', overwrite = TRUE)
terra::writeRaster(lc, 'Data/landCover.tiff', overwrite = TRUE)
terra::writeRaster(elevation, 'Data/elevation.tiff', overwrite = TRUE)

terra::writeRaster(Arable, 'Data/Arable.tiff', overwrite = TRUE)
terra::writeRaster(Coniferous, 'Data/Coniferous.tiff', overwrite = TRUE)
terra::writeRaster(Grassland, 'Data/Grassland.tiff', overwrite = TRUE)
terra::writeRaster(Built, 'Data/Built.tiff', overwrite = TRUE)

