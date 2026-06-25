ngroups <- 4
nEach <- 5
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

source('Misc/setupCovs.R')
Mesh <- readRDS('Data/MeshCoarse.rds')

BTOCV <- readRDS('Data/BTOCV.rds')

covsEx <- terra::extract(covariates_selection, BTOCV, ID = FALSE)

naRows <- lapply(covsEx, function(x) which(is.na(x)))
naCovs <- names(naRows)[sapply(naRows, length) > 0]

for(cov in naCovs){  # fill missing values for rows/covs using nearest neighbour
  covsEx[naRows[[cov]], cov] <-
    PointedSDMs::nearestValue(matrix(st_coordinates(BTOCV[naRows[[cov]],])[,c("X","Y")], ncol = 2),
                              terra::project(covariates_selection,
                                             crs)[cov])

}

precVals <- terra::extract(precVar, BTOCV, ID = FALSE, layer = as.character(BTOCV$year))
covsEx$Precipitation <- precVals$value
covsEx$year <- precVals$layer

for (year in min(BTOCV$year):2024) {

  precUsed <- subset(precVar, as.character(year))
  names(precUsed) <- 'Precipitation'

  naRows <- intersect(sapply(covsEx, function(x) which(is.na(x)))$Precipitation, which(covsEx$year == year))

  if (!identical(naRows, integer(0))){

  covsEx[naRows, 'Precipitation'] <-
    PointedSDMs::nearestValue(matrix(st_coordinates(BTOCV[naRows,])[,c("X","Y")], ncol = 2),
                              terra::project(precUsed,
                                             crs))

  }


}
covsEx$year <- NULL

BTOCV <- cbind(BTOCV, covsEx)

saveRDS(object = BTOCV, 'Data/BTOCV_Pred.rds')
