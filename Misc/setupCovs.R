crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

if (!exists('scaleVars')) scaleVars <- TRUE

library(INLA)
library(inlabru)
library(terra)
library(dplyr)
library(PointedSDMs)

if (!exists('timeDataTable')) {

  timeDataTable <- data.frame(Year = 2005:2024,
                            seq = rep(1:ngroups, each = nEach)[1:20])


#timeDataTable$seq <- c(rep(1, 5), rep(2, 6), rep(3, 6))

}
#timeDataTable$seq <- c(rep(1,8), rep(2, 9))

precVar <- terra::rast('Data/prec.tiff')
precVar2021 <- precVar[['2021']]
#precVar <- terra::rast('Data/precNew.tiff')
tempVar <- terra::rast('Data/temp.tiff')#tempBio4.tiff
elevVar <- terra::rast('Data/elevation.tiff')
#land <- terra::rast('Data/landCover.tiff')
Arable <- terra::rast('Data/Arable.tiff')
Coniferous <- terra::rast('Data/Coniferous.tiff')
Grassland <- terra::rast('Data/Grassland.tiff')
Built <- terra::rast('Data/Built.tiff')

BroadLeaf <- terra::rast('Data/Broadleaf.tiff')
SGrassland <- terra::rast('Data/SGrassland.tiff')
Mountain <- terra::rast('Data/Mountain.tiff')

Distance <- terra::rast('Data/Distance.tiff')


tempVar <- resample(tempVar, Arable)#land
elevVar <- resample(elevVar, Arable) #tempVar
precVar <- resample(precVar, Arable)
Arable <- resample(Arable, Arable)
Coniferous <- resample(Coniferous, Arable)
Grassland <- resample(Grassland, Arable)
Built <- resample(Built, Arable)


BroadLeaf <- resample(BroadLeaf, Arable)
SGrassland <- resample(SGrassland, Arable)
Mountain <- resample(Mountain, Arable)

Distance <- resample(Distance, Arable)

if (exists('timeDataTable')) {

  precVar <- terra::subset(precVar, subset = as.character(timeDataTable$Year))
  #tempVar <- terra::subset(tempVar, subset = as.character(timeDataTable$Year))


} else {

  precVar <- terra::subset(precVar, subset = as.character(2005:2024))
  tempVar <- terra::subset(tempVar, subset = as.character(2005:2021))

}

ElevVarSquared <- app(elevVar, function(x) x^2)
names(ElevVarSquared) <- 'elevationS'

precVarSquared <- app(precVar, function(x) x^2)

if (scaleVars) {

tempVals <- values(tempVar)
mean_valTemp <- mean(tempVals, na.rm = TRUE)
sd_valTemp   <- sd(tempVals, na.rm = TRUE)
tempVar <- (tempVar - mean_valTemp) / sd_valTemp

precVals <- values(precVar)
mean_valPrec <- mean(precVals, na.rm = TRUE)
sd_valPrec   <- sd(precVals, na.rm = TRUE)
precVar <- (precVar - mean_valPrec) / sd_valPrec

precValsS <- values(precVarSquared)
mean_valPrecS <- mean(precValsS, na.rm = TRUE)
sd_valPrecS   <- sd(precValsS, na.rm = TRUE)
precVarSquared <- (precVarSquared - mean_valPrecS) / sd_valPrecS

}

covariates_selection <- c(elevVar, ElevVarSquared, Arable, Coniferous, Grassland, Built,
                          BroadLeaf, SGrassland, Mountain, Distance)
if (scaleVars) covariates_selection <- scale(covariates_selection)
covariates_selection <- c(precVar$`2009`, tempVar$`2009`, precVarSquared$`2009`,covariates_selection)
#covariates_selection <- c(covariates_selection, land$`lc2018`)
#names(covariates_selection) <- c('Precipitation', 'Temperature', 'Elevation', 'LandCover')
names(covariates_selection) <- c('Precipitation', 'Temperature', 'PrecipitationS', 'Elevation', 'ElevationS','Arable', 'Coniferous', 'Grassland', 'Built',
                                 'BroadLeaf', 'SGrassland', 'Mountain', 'Distance')
#covariates_selection <- na.omit(covariates_selection)
covsUsed <- names(covariates_selection)
#covariates_selection <- covariates_selection[[1:3]]
covariates_selection <- terra::subset(covariates_selection, c('Temperature', 'Distance'), negate = TRUE) #Precipitation

UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)
UKsimp <- st_buffer(UKsimp, 5)
