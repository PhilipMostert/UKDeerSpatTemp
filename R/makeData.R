library(INLA)
library(inlabru)
library(terra)
library(dplyr)
library(PointedSDMs)

ngroups <- 4#8
nEach <-  4#2
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
species <- c("Capreolus capreolus", "Muntiacus reevesi", "Cervus elaphus", "Dama dama")
#Try Mesh coarse
speciesSelect <-  'Dama dama'
Mesh <- readRDS('Data/MeshFine.rds') #MeshFine.rds
timeDataTable <- data.frame(Year = 2005:2020,
                            seq = rep(1:ngroups, each = nEach))

# timeDataTable <- data.frame(Year = 2005:2020,
#                             seq = rep(1:2, each = 8))

precVar <- terra::rast('Data/prec.tiff')
tempVar <- terra::rast('Data/temp.tiff')#tempBio4.tiff
elevVar <- terra::rast('Data/elevation.tiff')
#land <- terra::rast('Data/landCover.tiff')
Arable <- terra::rast('Data/Arable.tiff')
Coniferous <- terra::rast('Data/Coniferous.tiff')
Grassland <- terra::rast('Data/Grassland.tiff')
Built <- terra::rast('Data/Built.tiff')

tempVar <- resample(tempVar, Arable)#land
elevVar <- resample(elevVar, Arable) #tempVar
precVar <- resample(precVar, Arable)
Arable <- resample(Arable, Arable)
Coniferous <- resample(Coniferous, Arable)
Grassland <- resample(Grassland, Arable)
Built <- resample(Built, Arable)

covariates_selection <- c(precVar$`2014`, tempVar$`2014`, elevVar, Arable, Coniferous, Grassland, Built)
covariates_selection <- scale(covariates_selection)
#covariates_selection <- c(covariates_selection, land$`lc2018`)
#names(covariates_selection) <- c('Precipitation', 'Temperature', 'Elevation', 'LandCover')
names(covariates_selection) <- c('Precipitation', 'Temperature', 'Elevation', 'Arable', 'Coniferous', 'Grassland', 'Built')
#covariates_selection <- na.omit(covariates_selection)
covsUsed <- names(covariates_selection)
#covariates_selection <- covariates_selection[[1:3]]
covariates_selection <- covariates_selection[[-1]] #No Precipitation
UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)
UKsimp <- st_buffer(UKsimp, 5)

biasField <- TRUE
tempBias <- TRUE
SVC <- FALSE
covariates_selection[['Elevation']] <- NULL
#covariates_selection[['Built']] <- NULL
#if (biasField) covariates_selection <- covariates_selection[[-length(names(covariates_selection))]]

for(spec in species) {

  speciesSelect <- spec

for (ind in 1:ngroups) {

print(ind)
BTO <- readRDS('data/BTO.rds')
iNat <- readRDS('data/iNat.rds')
Mammal <- readRDS('data/Mammal.rds')
MammalWeb <- readRDS('data/MammalWeb.rds')

PODat <- rbind(iNat,Mammal, MammalWeb)

PODat <- PODat %>% filter(species == speciesSelect, year %in% timeDataTable[timeDataTable$seq == ind,]$Year) %>% group_by(year) %>% distinct()
BTO <- BTO %>% filter(species == speciesSelect, year %in% timeDataTable[timeDataTable$seq == ind,]$Year) %>% mutate(listLength = scale(listLength))

#Ipoints <- fm_int(Mesh$Mesh, Mesh$UK, int.args = list(method = 'direct', nsub2 = 1, nsub1 = 1))
Ipoints <- fm_int(Mesh$Mesh, Mesh$UK)
Ipoints <- Ipoints[sapply(st_intersects(Ipoints,Mesh$UK), function(z) if (length(z)==0) FALSE else TRUE),]

model_setup <- startISDM(BTO, PODat, #PODat
                         Boundary = st_as_sf(Mesh$UK),
                         IPS = Ipoints,
                         Mesh = Mesh$Mesh,
                         Projection = crs,
                         spatialCovariates = covariates_selection,
                         responsePA = 'pres',
                         temporalName = 'year')#,
                         #pointCovariates = 'listLength'
                         #)

temp <- scale(tempVar)
prec <- scale(precVar)

prec <- prec[[as.character(timeDataTable[timeDataTable$seq == ind,]$Year)]]
temp <- temp[[as.character(timeDataTable[timeDataTable$seq == ind,]$Year)]]

#omit <- (filterYear1:filterYear2)[!filterYear1:filterYear2 %in% unique(PODat$year)]

for (i in seq(1, length(names(prec)))) {

  layer <- names(prec)[i]

  # model_setup$.__enclos_env__$private$IPS$Precipitation[model_setup$.__enclos_env__$private$IPS$year == i] <- eval_spatial(data = prec,
  #                                                                                                                           where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
  #                                                                                                                           layer = layer)
  #
  # model_setup$.__enclos_env__$private$IPS$Precipitation[ model_setup$.__enclos_env__$private$IPS$year == i] <- bru_fill_missing(data = prec, values =   model_setup$.__enclos_env__$private$IPS$Precipitation[ model_setup$.__enclos_env__$private$IPS$year == i],
  #                                                                                                                           where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
  #                                                                                                                           layer = layer)

  model_setup$.__enclos_env__$private$IPS$Temperature[ model_setup$.__enclos_env__$private$IPS$year == i] <- eval_spatial(data = temp,
                                                                                                                            where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
                                                                                                                            layer = layer)

  model_setup$.__enclos_env__$private$IPS$Temperature[ model_setup$.__enclos_env__$private$IPS$year == i] <- bru_fill_missing(data = temp, values =   model_setup$.__enclos_env__$private$IPS$Temperature[ model_setup$.__enclos_env__$private$IPS$year == i],
                                                                                                                                where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
                                                                                                                                layer = layer)

if (!is.null(model_setup$.__enclos_env__$private$modelData$BTO)) {



  model_setup$.__enclos_env__$private$modelData$BTO$BTO$Temperature[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i] <- eval_spatial(data = temp,
                                                                                                                          where = model_setup$.__enclos_env__$private$modelData$BTO$BTO[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i,],
                                                                                                                          layer = layer)

  model_setup$.__enclos_env__$private$modelData$BTO$BTO$Temperature[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i] <- bru_fill_missing(data = temp, values = model_setup$.__enclos_env__$private$modelData$BTO$BTO$Temperature[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i],
                                                                                                                              where = model_setup$.__enclos_env__$private$modelData$BTO$BTO[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i,],
                                                                                                                              layer = layer)

}

  if (!is.null(model_setup$.__enclos_env__$private$modelData$PODat)) {

  model_setup$.__enclos_env__$private$modelData$PODat$PODat$Temperature[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i] <- eval_spatial(data = temp,
                                                                                                                                                     where = model_setup$.__enclos_env__$private$modelData$PODat$PODat[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i,],
                                                                                                                                                     layer = layer)

  model_setup$.__enclos_env__$private$modelData$PODat$PODat$Temperature[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i] <- bru_fill_missing(data = temp, values = model_setup$.__enclos_env__$private$modelData$PODat$PODat$Temperature[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i],
                                                                                                                                                         where = model_setup$.__enclos_env__$private$modelData$PODat$PODat[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i,],
                                                                                                                                                         layer = layer)
}

}


if (biasField) {

if (!is.null(model_setup$.__enclos_env__$private$modelData$PODat)) {

  model_setup$addBias(datasetNames = 'PODat')

    if (!tempBias) model_setup$changeComponents('PODat_biasField(main = geometry, model = PODat_bias_field)')
    else {

      model_setup$changeComponents('PODat_biasField(main = geometry, model = PODat_bias_field, group = year, ngroup = 4, control.group = list(model = "iid", hyper = list(prec = list(fixed = TRUE))))')
      #model_setup$changeComponents('PODat_spatial(main = geometry, copy = "BTO_spatial", group = year, control.group = list(model = "ar1"), hyper = list(beta = list(fixed = TRUE)))')
      #Try making temporal Bias iid

    }
  }

}


if (SVC) model_setup$changeComponents('Temperature(main = geometry, weights = Temperature, model = BTO_field)')

#model_setup$updateFormula(datasetName = 'BTO', Formula = ~ . - Built)


##Change these for Dama Dama and for bias?
#Lower kappa for smoother field
#Increase tau = rougher field

rho0 <- Mesh$Mesh$loc[c(Mesh$Mesh$segm$bnd$idx[,1], Mesh$Mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5
sigma0 <- 1
alpha <- 2; d <- 2; nu <- alpha - d/2
kappa0 <- log(8*nu)/2 - log(rho0)
tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0
spde <- inla.spde2.matern(mesh = Mesh$Mesh, B.tau = cbind(tau0, nu, -1),
                          #theta.prior.mean = c(1.8, 1), theta.prior.prec = c(10, 1), #1.8, 1
                          B.kappa = cbind(kappa0, -1, 0), alpha = alpha, constr = FALSE)

#
if (!is.null(model_setup$.__enclos_env__$private$modelData$PODat)) model_setup$spatialFields$datasetFields$BTO <- spde
else model_setup$spatialFields$sharedField$sharedField <- spde
if (biasField) model_setup$spatialFields$biasFields$PODat <- spde

print("Model setup done")
speciesPrint <- gsub(' ', '', speciesSelect)
saveRDS(model_setup, paste0(speciesPrint,'/', speciesPrint,'Sub', ind,'.rds'))

}
print(spec)

}

