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

for (i in 1:length(list.files('Climate/tempNew/'))) {

file <- list.files(path = 'Climate/tempNew/')[i]
xx <- terra::rast(paste0('Climate/tempNew/', file))

xx <- mask(crop(xx, st_transform(UK, fm_wkt(xx))), st_transform(UK, fm_wkt(xx)))

if (i == 1) temp <- xx
else temp <- c(temp, xx)

}

names(temp) <- 2000:2024
temp <- project(temp, crs)

for (i in 1:length(list.files('Climate/precNew/'))) {

  file <- list.files(path = 'Climate/precNew/')[i]
  xx <- terra::rast(paste0('Climate/precNew/', file))

  xx <- mask(crop(xx, st_transform(UK, fm_wkt(xx))), st_transform(UK, fm_wkt(xx)))

  if (i == 1) prec <- xx
  else prec <- c(prec, xx)

}

names(prec) <- 2000:2024
prec <- project(prec, crs)

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
BroadLeaf <- project(LC$Layer_1, crs(UK))
SGrassland <- project(LC$Layer_5, crs(UK))
Mountain <- project(LC$Layer_6, crs(UK))


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

BroadLeaf <- mask(crop(BroadLeaf, UK), UK)
BroadLeaf <- project(BroadLeaf, crs)
names(BroadLeaf) <- 'BroadLeaf'

SGrassland <- mask(crop(SGrassland, UK), UK)
SGrassland <- project(SGrassland, crs)
names(SGrassland) <- 'SGrassland'

Mountain <- mask(crop(Mountain, UK), UK)
Mountain <- project(Mountain, crs)
names(Mountain) <- 'Mountain'


covs <- list(Temperature = temp,
             Precipitation = prec,
             LandCover = lc,
             Elevation = elevation)

terra::writeRaster(temp, 'Data/temp.tiff', overwrite = TRUE)
terra::writeRaster(prec, 'Data/prec.tiff', overwrite = TRUE)
terra::writeRaster(elevation, 'Data/elevation.tiff', overwrite = TRUE)

terra::writeRaster(Arable, 'Data/Arable.tiff', overwrite = TRUE)
terra::writeRaster(Coniferous, 'Data/Coniferous.tiff', overwrite = TRUE)
terra::writeRaster(Grassland, 'Data/Grassland.tiff', overwrite = TRUE)
terra::writeRaster(Built, 'Data/Built.tiff', overwrite = TRUE)

terra::writeRaster(BroadLeaf, 'Data/BroadLeaf.tiff', overwrite = TRUE)
terra::writeRaster(SGrassland, 'Data/SGrassland.tiff', overwrite = TRUE)
terra::writeRaster(Mountain, 'Data/Mountain.tiff', overwrite = TRUE)

