print('Starting')
library(INLA)
library(inlabru)
library(terra)
library(dplyr)
library(PointedSDMs)

timeParIn <-  20
yearsIn <- 2005:2024

args <- commandArgs(trailingOnly = TRUE)

File <- args[1]

results <- readRDS(paste0(File, 'modEsts.rds'))

latent_fixed_effects <- results$latent_fixed_effects
x <- results$x
rm(results)

numTimes <- length(latent_fixed_effects)

print(numTimes)
print(File)
model <- readRDS(paste0(File, 'finalModEsts.rds'))

thetaMode <- model$thetaMode

intDesign <- model$intDesign

modNames <- row.names(latent_fixed_effects[[1]])

spdeMain <- model$spdeMain
spdeBias <- model$spdeBias
spdeCommon <- model$spdeCommon
rm(model)

print('Estimating')

predInd <- rep(1:numTimes, each = timeParIn/numTimes)
BTOCV <- readRDS('Data/BTOCV_Pred.rds')

SpecUsed <- switch(gsub("Sub1.rds$", "", sort(list.files(paste0(File, 'Data/')))[1]),
                   Muntiacusreevesi = 'Muntiacus reevesi',
                   Damadama = 'Dama dama',
                   Capreoluscapreolus = 'Capreolus capreolus',
                   Cervuselaphus = 'Cervus elaphus')

print(SpecUsed)

BTO_Pred <- BTOCV[BTOCV$species == SpecUsed, ]

CVListLog <- CVList <- list()
DIC <- list()
WAIC <- list()

if (!dir.exists(paste0(File, 'Summaries'))) dir.create(paste0(File, 'Summaries'))

for(i in seq(1, numTimes)) {
print(i)

file <- grep(i, list.files(paste0(File, '/Data/')))
file <- list.files(paste0(File, '/Data/'))[file][1]
model_setup <- readRDS(paste0(File, '/Data/', file))

print(file)

BTO_Temp <- BTO_Pred[BTO_Pred$year %in% (yearsIn)[predInd == i], ]
BTO_Temp$yearTemp <- BTO_Temp$year
BTO_Temp$year <- as.numeric(as.factor(rank(BTO_Temp$year)))


if(i==1){

  fixed.mean_corrected <- ((latent_fixed_effects[[numTimes]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2)**(-1)*(latent_fixed_effects[[numTimes]]$sd**(-2)*latent_fixed_effects[[numTimes]]$mean - latent_fixed_effects[[i]]$sd**(-2)*latent_fixed_effects[[i]]$mean)) #%>% list(.)
  fixed.prec_corrected <- (latent_fixed_effects[[numTimes]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2)# %>% list(.)
  names(fixed.mean_corrected) <- names(fixed.prec_corrected) <- modNames#<- rownames(model$summary.fixed)

} else if(i== numTimes){

  fixed.mean_corrected <- latent_fixed_effects[[i-1]]$mean ##%>% list(.)
  fixed.prec_corrected <- latent_fixed_effects[[i-1]]$sd**-2 ##%>% list(.)
  names(fixed.mean_corrected) <- names(fixed.prec_corrected) <- modNames#<- rownames(model$summary.fixed)

} else{

  fixed.mean_corrected <- ((latent_fixed_effects[[i-1]]$sd**-2 + latent_fixed_effects[[numTimes]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2)**(-1)*(latent_fixed_effects[[i-1]]$sd**(-2)*latent_fixed_effects[[i-1]]$mean + latent_fixed_effects[[numTimes]]$sd**(-2)*latent_fixed_effects[[numTimes]]$mean - latent_fixed_effects[[i]]$sd**(-2)*latent_fixed_effects[[i]]$mean)) #%>% list(.)
  fixed.prec_corrected <- (latent_fixed_effects[[i-1]]$sd**-2 + latent_fixed_effects[[numTimes]]$sd**-2 - latent_fixed_effects[[i]]$sd**-2)# %>% list(.)
  names(fixed.mean_corrected) <- names(fixed.prec_corrected) <- modNames#<- rownames(model$summary.fixed)

}

  for (cov in names(fixed.mean_corrected)[!grepl('intercept', names(fixed.mean_corrected))]) {
    model_setup$priorsFixed(cov, mean.linear = fixed.mean_corrected[cov],
                            prec.linear = fixed.prec_corrected[cov])
  }

if (all(c('BTO', 'PODat') %in% model_setup$.__enclos_env__$private$dataSource)) {

  for (int in c('BTO', 'PODat')) {

    model_setup$priorsFixed(Effect = 'Intercept',
                            mean.linear = fixed.mean_corrected[paste0(int,'_intercept')],
                            prec.linear = fixed.prec_corrected[paste0(int, '_intercept')],
                            datasetName = int)
  }

   model_setup$spatialFields$datasetFields[[1]] <- spdeMain
   model_setup$spatialFields$biasFields$PODat <- spdeBias
   if (!is.null(model_setup$spatialFields$biasFields$sharedBias)) model_setup$spatialFields$biasFields$sharedBias <- spdeCommon

} else {

      model_setup$priorsFixed(Effect = 'Intercept',
                              mean.linear = fixed.mean_corrected[paste0('BTO','_intercept')],
                              prec.linear = fixed.prec_corrected[paste0('BTO', '_intercept')],
                              datasetName = 'BTO')

      model_setup$spatialFields$sharedField$sharedField <- spdeMain
      if (!is.null(model_setup$spatialFields$biasFields$sharedBias)) model_setup$spatialFields$biasFields$sharedBias <- spdeCommon

}


 YEARIND <- grepl('yearTrend\\(main = year, model = "ar1"',model_setup$.__enclos_env__$private$Components)

 if (any(YEARIND)) {

 if (i > 1) {

       model_setup$.__enclos_env__$private$Components[YEARIND] <- 'yearTrend(main = year, model = "ar1", constr = FALSE,hyper = list(mean = list(initial = MODEVAR), prec = list(prior = "gaussian", param = c(0, 1))))'

 }

}

  model_setup$.__enclos_env__$private$Components <- sort(model_setup$.__enclos_env__$private$Components)


  model_rf <- fitISDM(model_setup, options = list(num.threads = 4,
                                                  safe = FALSE,
                                                  bru_max_iter = 1,
                                                  inla.mode = 'experimental',
                                                  control.mode = list(theta = thetaMode),
                                                  control.compute = list(dic = TRUE, waic = TRUE),
                                                  control.inla = list(
                                                    cmin = 0,
                                                    control.vb = list(enable = FALSE),
                                                    int.strategy = 'user.expert',
                                                    int.design = intDesign,
                                                    strategy = "adaptive",
                                                    diagonal = 0)))


  DIC[[i]] <- model_rf$dic
  WAIC[[i]] <- model_rf$waic
  POPA <- all(c('PODat', 'BTO') %in% model_setup$.__enclos_env__$private$dataSource)
  whichVar <-  c('year', 'yearTemp', 'pres', 'pres2', 'n_visits', 'species','mean', 'sd',  'q0.025', 'q0.5', 'q0.975')

  fixedComps <- row.names(model_rf$summary.fixed)[!row.names(model_rf$summary.fixed) %in% c('Distance', 'PODat_intercept', 'transects_year')]
  randomComps <- names(model_rf$summary.random)[!names(model_rf$summary.random) %in% c('POTrend', 'PODat_spatial', 'PODat_biasField')]

  predictions <- data.frame(generate(object = model_rf, newdata = BTO_Temp, formula = formula(paste('~', paste0(c(fixedComps, randomComps), collapse = '+'))), n.samples = 1000))

  colnames(predictions) <- paste0('Sample:', 1:ncol(predictions))
  predictions$year <- BTO_Temp$yearTemp
  predictions$pres <- BTO_Temp$pres

  CVListLog[[i]] <- predictions


  for (j in 1:length(model_rf$misc$configs$config)) {

        model_rf$misc$configs$config[[j]]$Qinv <- NULL
        model_rf$misc$config$config[[j]]$Qprior <- NULL
        model_rf$misc$config$config[[j]]$mean <- NULL
        model_rf$misc$config$config[[j]]$ll.info <- NULL
        model_rf$misc$config$config[[j]]$Predictor <- NULL

  }

  model_rf$.args = NULL
  model_rf$summary.fitted.values <- NULL
  model_rf$.args <- NULL
  model_rf$summary.linear.predictor <- NULL
  model_rf$optionsJoint <- NULL
  model_rf$bru_iinla <- NULL
  model_rf$marginals.random <- NULL


saveRDS(model_rf$summary.random, paste0(File,'Summaries/sumRandom', i,'.rds'))
saveRDS(model_rf$summary.hyper, paste0(File,'Summaries/sumHyper', i, '.rds'))
saveRDS(model_rf$summary.fixed, paste0(File,'Summaries/fixedRes', i, '.rds'))


if (i == numTimes) {


if (!is.null(model_rf$bru_info$model$effects[['PODat_spatial']]$main$model)) saveRDS(INLA::inla.spde2.result(inla = model_rf,
                                name = 'PODat_spatial',
                                spde = model_rf$bru_info$model$effects[['PODat_spatial']]$main$model),
        file = paste0(File, 'Summaries/PODat_spatial.rds'))
 else saveRDS(INLA::inla.spde2.result(inla = model_rf,
                                      name = 'shared_spatial',
                                      spde = model_rf$bru_info$model$effects[['shared_spatial']]$main$model),
              file = paste0(File, 'Summaries/BTOResult.rds'))


if (!is.null(model_rf$bru_info$model$effects[['PODat_biasField']]$main$model)) saveRDS(INLA::inla.spde2.result(inla = model_rf,
                                name = 'PODat_biasField',
                                spde = model_rf$bru_info$model$effects[['PODat_biasField']]$main$model),
        file = paste0(File, 'Summaries/BiasResult.rds'))

 marginal_trend <- model_rf$internal.marginals.hyperpar$`Log precision for yearTrend`

 saveRDS(marginal_trend, paste0(File, 'Summaries/trendResult.rds'))

 if (!is.null(model_rf$internal.marginals.hyperpar$`Log precision for POTrend`)) {

  marginal_POtrend <- model_rf$internal.marginals.hyperpar$`Log precision for POTrend`
  saveRDS(marginal_POtrend, paste0(File, 'Summaries/POTrendResult.rds'))


 }

}

model_rf$misc <- NULL

MODEVAR <- tail(model_rf$summary.random$yearTrend$mode, 1)
MODEVARPO <- tail(model_rf$summary.random$POTrend$mode, 1)

print(paste("Model", i, "is done."))
rm(model_rf)

}

CVLog <- do.call(rbind, CVListLog)

saveRDS(DIC, paste0(File, 'DIC.rds'))
saveRDS(WAIC, paste0(File, 'WAIC.rds'))
saveRDS(CVLog, paste0(File, 'intensitySamps.rds'))

print(File)
