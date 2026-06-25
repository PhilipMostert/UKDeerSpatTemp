library(terra)
library(sf)
library(inlabru)
library(INLA)
library(ggplot2)
library(tidyterra)

save <- TRUE
yearsIn <- 2005:2024

Folder <- 'xx/'

speciesIn <- c('Capreoluscapreolus', 'Cervuselaphus', 'Damadama', 'Muntiacusreevesi')

Mesh <- readRDS('Data/MeshCoarse.rds')

for (speciesName in speciesIn) {
print(speciesName)


nIn <- length(list.files(paste0(Folder, speciesName, '/modelRes/')))

if (!dir.exists(paste0(Folder, speciesName, '/Results')))  dir.create(paste0(Folder, speciesName, '/Results'))

predData <- readRDS('Data/predData1kmList.rds')

coords <- data.frame(st_coordinates(predData[[1]]))
names(coords) <- c('x', 'y')

Random <- list()

for (i in 1:nIn) {

  RFile <- readRDS(paste0(Folder, speciesName,'/Summaries/sumRandom', i,'.rds'))

  if (i == 1) {

    Random$PODat_biasField <- RFile$PODat_biasField

    if (!is.null(RFile$PODat_spatial)) Random$BTO_spatial <- RFile$BTO_spatial
    else Random$BTO_spatial <- RFile$shared_spatial

    Random$yearTrend <- RFile$yearTrend
    Random$POTrend <- RFile$POTrend

  } else {

    Random$PODat_biasField <- rbind(Random$PODat_biasField, RFile$PODat_biasField)

    if (!is.null(RFile$PODat_spatial))  Random$BTO_spatial <- rbind(Random$BTO_spatial, RFile$BTO_spatial)
    else Random$BTO_spatial <- rbind(Random$BTO_spatial, RFile$shared_spatial)

    Random$yearTrend <- rbind(Random$yearTrend, RFile$yearTrend)
    Random$POTrend <- rbind(Random$POTrend, RFile$POTrend)


  }

}


Fixed <- readRDS(paste0(Folder, speciesName,'/Summaries/fixedRes', nIn,'.rds'))

FixedList <- list()

FixedUsed <- Fixed[row.names(Fixed) %in% names(predData[[1]]),]

for (i in 1:length(predData)) {

  Vals <- st_drop_geometry(predData[[i]])
  Vals <- Vals[,names(Vals) %in% row.names(Fixed)]

  Vals <- data.frame(as.matrix(Vals) %*% as.matrix(FixedUsed[names(Vals),]))

  linPred <- cbind(coords, data.frame(median = Vals$X0.5quant))
  linPred <- rast(linPred)

crs(linPred) <- fm_wkt(Mesh$Mesh)
linPred <- mask(linPred, Mesh$UK)
FixedList[[i]] <- linPred
names(FixedList)[[i]] <- as.character(i)

}

if (!is.null(Random$PODat_biasField)) {

if (nrow(Random$PODat_biasField) == Mesh$Mesh$n * length(yearsIn)) {

BiasDat <- split(Random$PODat_biasField, cut(seq_along(Random$PODat_biasField$mean), length(yearsIn), labels = FALSE))

mean_bias <- rast()
sd_bias <- rast()

for (i in seq(1, length(yearsIn))) {

  BiasField <- fm_evaluate(Mesh$Mesh, data.frame(BiasDat[[i]]), loc = predData[[1]])
  mean_j <- cbind(coords, data.frame(BiasField[, c('sd', 'X0.5quant')]))

  names(mean_j)[names(mean_j) == 'X0.5quant'] <- 'median'

  mean_j <- rast(mean_j)
  crs(mean_j) <- fm_wkt(Mesh$Mesh)
  mean_j <- mask(mean_j, Mesh$UK)

  median <- mean_j$median
  names(median) <- yearsIn[i]
  sd <- mean_j$sd
  names(sd) <- yearsIn[i]

  mean_bias <- c(mean_bias, median)
  sd_bias <- c(sd_bias, sd)
}

ggplot() + tidyterra::geom_spatraster(data = mean_bias) +
  facet_wrap(~lyr, nrow = 2) +
  scale_fill_grass_c(palette = 'inferno',  use_grass_range = FALSE) +
  theme_minimal()

  if (save) {

    saveRDS(object = mean_bias, file = paste0(Folder, speciesName ,'/Results/', speciesName,'_biasMedian', '.rds'))
    saveRDS(object = sd_bias, file = paste0(Folder, speciesName ,'/Results/', speciesName,'_biasSD', '.rds'))


  }




} else {

  BiasField <- fm_evaluate(Mesh$Mesh, RFile$PODat_biasField, loc = predData[[1]])
  mean_j <- cbind(coords, data.frame(BiasField[, c('sd', '0.5quant')]))
  names(mean_j)[names(mean_j) == 'X0.5quant'] <- 'median'
  mean_j <- rast(mean_j)
  crs(mean_j) <- fm_wkt(Mesh$Mesh)
  mean_bias <- mask(mean_j, Mesh$UK)
  plot(mean_bias)
  if(save) saveRDS(object = mean_bias, file = paste0(Folder, speciesName ,'/Results/', speciesName,'_biasPreds', '.rds'))

}
}

TimeDat <- split(Random$BTO_spatial, cut(seq_along(Random$BTO_spatial$mean), length(yearsIn), labels = FALSE))
mean_all <- rast()
sd_all <- rast()

for (i in seq(1, length(yearsIn))) {

  TimeField <- fm_evaluate(Mesh$Mesh,  data.frame(TimeDat[[i]]), loc = predData[[i]])
  mean_j <- cbind(coords,  data.frame(TimeField[, c('sd', 'X0.5quant')]))
  names(mean_j)[names(mean_j) == 'X0.5quant'] <- 'median'
  mean_j <- rast(mean_j)
  crs(mean_j) <- fm_wkt(Mesh$Mesh)
  mean_j <- mask(mean_j, Mesh$UK)
  median <- mean_j$median + Random$yearTrend[['0.5quant']][i] + FixedList[[i]]$median + Fixed['BTO_intercept', '0.5quant']
  sd <- mean_j$sd + Random$yearTrend[['sd']][i]

  names(median) <-  yearsIn[i]
  names(sd) <- yearsIn[i]

  mean_all <- c(mean_all, median)
  sd_all <- c(sd_all, sd)

}

if (save) {

  saveRDS(object = mean_all,  file = paste0(Folder, speciesName, '/Results/', speciesName, '_linPredMedian', '.rds'))
  saveRDS(object = sd_all,  file = paste0(Folder, speciesName, '/Results/', speciesName, '_linPredSD', '.rds'))

}


}
