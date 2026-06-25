library(PointedSDMs)
library(ggplot2)
library(terra)
library(inlabru)

ListVals <- TRUE

crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)
UKsimp <- st_buffer(UKsimp, 5)

yearsIn <- as.character(2005:2024)


Mesh <- readRDS('Data/MeshCoarse.rds')$Mesh

bbox <- fm_bbox(UKsimp)

dim <- 1

dimsOut <- round(c(x = diff(bbox[[1]]), y = diff(bbox[[2]])) / dim)

#Create predData from projectorMesh
projectorMesh <- fm_evaluator(Mesh,
                              xlim = bbox[[1]], ylim = bbox[[2]], dims = dimsOut)

predData <- st_as_sf(expand.grid(x = projectorMesh$x, y = projectorMesh$y),
                     coords = c('x', 'y'),
                     crs = fm_crs(UKsimp, crs))

whichIn <- unlist(st_contains(st_transform(UKsimp, crs), predData))
predData <- predData[whichIn, ]
ngroups <- 4
nEach <-  5
source('Misc/setupCovs.R')

covs <- terra::extract(covariates_selection, predData, ID = FALSE)

naRows <- lapply(covs, function(x) which(is.na(x)))
naCovs <- names(naRows)[sapply(naRows, length) > 0]

for(cov in naCovs){  # fill missing values for rows/covs using nearest neighbour
  covs[naRows[[cov]], cov] <-
    PointedSDMs::nearestValue(matrix(st_coordinates(predData[naRows[[cov]],])[,c("X","Y")], ncol = 2),
                 terra::project(covariates_selection,
                                crs)[cov])

}

precVals <- values(precVar)
mean_valPrec <- mean(precVals, na.rm = TRUE)
sd_valPrec   <- sd(precVals, na.rm = TRUE)
precVar <- (precVar - mean_valPrec) / sd_valPrec
precVar <- resample(precVar, Arable)

precList <- list()

if (ListVals) {

  covs$Precipitation <- covs$PrecipitationS <- NULL

for (year in 1:length(yearsIn)) {

  precUsed <- precVar[[year]]
  names(precUsed) <- 'Precipitation'

  covPrec <- terra::extract(precUsed, predData, ID = FALSE)

  naRows <- sapply(covPrec, function(x) which(is.na(x)))

  covPrec[naRows, ] <-
      PointedSDMs::nearestValue(matrix(st_coordinates(predData[naRows,])[,c("X","Y")], ncol = 2),
                                terra::project(precUsed,
                                               crs))

   precList[[yearsIn[year]]] <- cbind(predData, covs, covPrec)


}

  saveRDS(precList, file = 'Data/predData1kmList.rds')

} else {

predData <- cbind(predData, covs)
saveRDS(predData, 'Data/predData1km.rds')

}
