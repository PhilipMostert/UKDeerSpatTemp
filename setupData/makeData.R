library(PointedSDMs)
library(INLA)
biasField <- TRUE
tempBias <-  TRUE

tempTrend <- TRUE

rhoPrior <- TRUE

Quads <- FALSE
iidTemp <- FALSE

ngroups <- 4
nEach <-  5

timeDataTable <- data.frame(Year = 2005:2024, #2009
                            seq = rep(1:ngroups, each = nEach)[1:20], #16
                            ind = 1:20)#16

crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
species <- c( "Muntiacus reevesi", "Dama dama","Capreolus capreolus", "Cervus elaphus")

Mesh <- readRDS('Data/MeshCoarse.rds') #MeshFine.rds

source('Misc/setupCovs.R')


if (!Quads) covariates_selection <- terra::subset(covariates_selection, c('PrecipitationS', 'ElevationS'), negate = TRUE)

for(spec in species) {

  speciesSelect <- spec
  speciesPrint <- gsub(' ', '', speciesSelect)

  if (!dir.exists(paste0('SpeciesModelData/',speciesPrint, '/Data'))) dir.create(paste0('SpeciesModelData/',speciesPrint, '/Data'))

for (ind in 1:ngroups) {


print(ind)

BTO <- readRDS('data/BTO.rds')

iNat <- readRDS('data/iNat.rds')
Mammal <- readRDS('data/Mammal.rds')

PODat <- rbind(iNat, Mammal)

PODat <- PODat %>% filter(species == speciesSelect, year %in% timeDataTable[timeDataTable$seq == ind,]$Year) %>% group_by(year) %>% distinct()
BTO <- BTO %>% mutate(
                      transects_year = scale(log(transects_year))) %>% filter(species == speciesSelect, year %in% timeDataTable[timeDataTable$seq == ind,]$Year)

if ('IPS' %in% names(Mesh)) Ipoints <- Mesh$IPS
else {
if (Mesh$Mesh$n < 3000) Ipoints <- fm_int(Mesh$Mesh, Mesh$UK, int.args = list(method = 'direct', nsub2 = 1, nsub1 = 3))
else Ipoints <- fm_int(Mesh$Mesh, Mesh$UK)
}

model_setup <- startISDM(PODat,
                         BTO,
                         Boundary = st_as_sf(Mesh$UK),
                         IPS = Ipoints,
                         Mesh = Mesh$Mesh,
                         Projection = crs,
                         spatialCovariates = covariates_selection,
                         responsePA = 'pres',
                         pointCovariates = 'transects_year',
                         temporalName = 'year')


prec <- precVar[[as.character(timeDataTable[timeDataTable$seq == ind,]$Year)]]
precS <- precVarSquared[[as.character(timeDataTable[timeDataTable$seq == ind,]$Year)]]

if (any(c('Precipitation', 'Temperature') %in% names(covariates_selection))) {

for (i in seq(1, length(names(prec)))) {

  layer <- names(prec)[i]

  model_setup$.__enclos_env__$private$IPS$Precipitation[model_setup$.__enclos_env__$private$IPS$year == i] <- eval_spatial(data = prec,
                                                                                                                            where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
                                                                                                                            layer = layer)

  model_setup$.__enclos_env__$private$IPS$Precipitation[ model_setup$.__enclos_env__$private$IPS$year == i] <- bru_fill_missing(data = prec, values =   model_setup$.__enclos_env__$private$IPS$Precipitation[ model_setup$.__enclos_env__$private$IPS$year == i],
                                                                                                                            where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
                                                                                                                            layer = layer)

if ('PrecipitationS' %in% names(covariates_selection)) {

  model_setup$.__enclos_env__$private$IPS$PrecipitationS[model_setup$.__enclos_env__$private$IPS$year == i] <- eval_spatial(data = precS,
                                                                                                                           where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
                                                                                                                           layer = layer)

  model_setup$.__enclos_env__$private$IPS$PrecipitationS[ model_setup$.__enclos_env__$private$IPS$year == i] <- bru_fill_missing(data = precS, values =   model_setup$.__enclos_env__$private$IPS$PrecipitationS[ model_setup$.__enclos_env__$private$IPS$year == i],
                                                                                                                                where = model_setup$.__enclos_env__$private$IPS[ model_setup$.__enclos_env__$private$IPS$year == i,],
                                                                                                                                layer = layer)


}



if (!is.null(model_setup$.__enclos_env__$private$modelData$BTO)) {

  model_setup$.__enclos_env__$private$modelData$BTO$BTO$Precipitation[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i] <- eval_spatial(data = prec,
                                                                                                                                                     where = model_setup$.__enclos_env__$private$modelData$BTO$BTO[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i,],
                                                                                                                                                     layer = layer)

  model_setup$.__enclos_env__$private$modelData$BTO$BTO$Precipitation[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i] <- bru_fill_missing(data = prec, values = model_setup$.__enclos_env__$private$modelData$BTO$BTO$Precipitation[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i],
                                                                                                                                                         where = model_setup$.__enclos_env__$private$modelData$BTO$BTO[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i,],
                                                                                                                                                         layer = layer)

  if ('PrecipitationS' %in% names(covariates_selection)) {

    model_setup$.__enclos_env__$private$modelData$BTO$BTO$PrecipitationS[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i] <- eval_spatial(data = precS,
                                                                                                                                                         where = model_setup$.__enclos_env__$private$modelData$BTO$BTO[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i,],
                                                                                                                                                         layer = layer)

    model_setup$.__enclos_env__$private$modelData$BTO$BTO$PrecipitationS[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i] <- bru_fill_missing(data = precS, values = model_setup$.__enclos_env__$private$modelData$BTO$BTO$PrecipitationS[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i],
                                                                                                                                                             where = model_setup$.__enclos_env__$private$modelData$BTO$BTO[model_setup$.__enclos_env__$private$modelData$BTO$BTO$year == i,],
                                                                                                                                                             layer = layer)


    }



}

  if (!is.null(model_setup$.__enclos_env__$private$modelData$PODat)) {


    model_setup$.__enclos_env__$private$modelData$PODat$PODat$Precipitation[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i] <- eval_spatial(data = prec,
                                                                                                                                                               where = model_setup$.__enclos_env__$private$modelData$PODat$PODat[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i,],
                                                                                                                                                               layer = layer)

    model_setup$.__enclos_env__$private$modelData$PODat$PODat$Precipitation[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i] <- bru_fill_missing(data = prec, values = model_setup$.__enclos_env__$private$modelData$PODat$PODat$Precipitation[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i],
                                                                                                                                                                   where = model_setup$.__enclos_env__$private$modelData$PODat$PODat[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i,],
                                                                                                                                                                   layer = layer)

    if ('PrecipitationS' %in% names(covariates_selection)) {

      model_setup$.__enclos_env__$private$modelData$PODat$PODat$PrecipitationS[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i] <- eval_spatial(data = precS,
                                                                                                                                                                   where = model_setup$.__enclos_env__$private$modelData$PODat$PODat[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i,],
                                                                                                                                                                   layer = layer)

      model_setup$.__enclos_env__$private$modelData$PODat$PODat$PrecipitationS[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i] <- bru_fill_missing(data = precS, values = model_setup$.__enclos_env__$private$modelData$PODat$PODat$PrecipitationS[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i],
                                                                                                                                                                       where = model_setup$.__enclos_env__$private$modelData$PODat$PODat[model_setup$.__enclos_env__$private$modelData$PODat$PODat$year == i,],
                                                                                                                                                                       layer = layer)


    }

}

}

}

if (biasField) {

if (!is.null(model_setup$.__enclos_env__$private$modelData$PODat)) {

  model_setup$addBias(datasetNames = 'PODat')

    if (!tempBias) model_setup$changeComponents('PODat_biasField(main = geometry, model = PODat_bias_field)')
    else model_setup$changeComponents('PODat_biasField(main = geometry, model = PODat_bias_field, replicate = year)')

  }

}


if (tempTrend) {

  model_setup$changeComponents('yearTrend(main = year, model = "ar1", constr = FALSE, hyper = list(prec = list(prior = "gaussian", param = c(0, 2)), rho = list(prior ="gaussian", param = c(0,0.15))))')

    model_setup$updateFormula(datasetName = 'BTO',Formula = ~ . + yearTrend)
    if ('PODat' %in% model_setup$.__enclos_env__$private$dataSource) model_setup$updateFormula(datasetName = 'PODat',Formula = ~ . + yearTrend)

}

if (rhoPrior) {

if ('PODat' %in% model_setup$.__enclos_env__$private$dataSource) model_setup$changeComponents(paste0('PODat_spatial(main = geometry, model = PODat_field, group = year,
    ngroup = ', length(timeDataTable[timeDataTable$seq == ind,]$Year),', control.group = list(model = "ar1", hyper = list(theta = list(prior = "gaussian", param = c(0,0.5)))))'))
  else model_setup$changeComponents(paste0('shared_spatial(main = geometry, model = shared_field, group = year,
    ngroup = ', length(timeDataTable[timeDataTable$seq == ind,]$Year),', control.group = list(model = "ar1", hyper = list(theta = list(prior = "gaussian", param = c(0,0.5)))))'))

}

if (iidTemp) {

  if ('PODat' %in% model_setup$.__enclos_env__$private$dataSource) {

  model_setup$changeComponents(paste0('PODat_spatial(main = geometry, model = PODat_field, group = year,
    ngroup = ', length(timeDataTable[timeDataTable$seq == ind,]$Year),', control.group = list(model = "iid"))'))

    model_setup$changeComponents(paste0('BTO_spatial(main = geometry, copy = "PODat_spatial", group = year,
    control.group = list(model = "iid"), hyper = list(beta = list(fixed = TRUE)))'))

  }
  else model_setup$changeComponents(paste0('shared_spatial(main = geometry, model = shared_field, group = year,
    ngroup = ', length(timeDataTable[timeDataTable$seq == ind,]$Year),', control.group = list(model = "iid"))'))

}


##Add priors
rho0 <- Mesh$Mesh$loc[c(Mesh$Mesh$segm$bnd$idx[,1], Mesh$Mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5
sigma0 <- 1
alpha <- 2; d <- 2; nu <- alpha - d/2
kappa0 <- log(8*nu)/2 - log(rho0)
tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0


spdeSpatTemp <- inla.spde2.matern(mesh = Mesh$Mesh, B.tau = cbind(tau0, nu, -1),
                          B.kappa = cbind(kappa0, -1, 0),
                          alpha = alpha,
                          prior.variance.nominal = 0.1,
                          constr = TRUE)

if (!is.null(model_setup$.__enclos_env__$private$modelData$PODat)) model_setup$spatialFields$datasetFields[[1]] <- spdeSpatTemp
else model_setup$spatialFields$sharedField$sharedField <- spdeSpatTemp


if (biasField) {

  alpha <- 2; d <- 2; nu <- alpha - d/2
  kappa0 <- log(8*nu)/2 - log(rho0)
  tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0

spdeBias <- inla.spde2.matern(mesh = Mesh$Mesh, B.tau = cbind(tau0, nu, -1),
                              B.kappa = cbind(kappa0, -1, 0), alpha = alpha,
                              prior.variance.nominal = 0.1,
                              constr = ifelse(tempBias, FALSE, TRUE))

model_setup$spatialFields$biasFields$PODat <- spdeBias

}

print("Model setup done")

saveRDS(model_setup, paste0('SpeciesModelData/',speciesPrint,'/Data/', speciesPrint,'Sub', ind,'.rds'))

}
print(spec)

}
