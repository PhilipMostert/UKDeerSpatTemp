library(INLA)
library(inlabru)
library(terra)
library(dplyr)
library(PointedSDMs)

crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
species <- c("Capreolus capreolus", "Muntiacus reevesi", "Cervus elaphus", "Dama dama")

args <- commandArgs(trailingOnly = TRUE)

File <- args[1]

Obj1 <- list.files(paste0(File, 'Data/'))[grepl('1', list.files(paste0(File, 'Data/')))]
print(Obj1)
model_setup <- readRDS(paste0(File, 'Data/', Obj1[1]))

ar1Mod <- any(grepl('yearTrend\\(main = year, model = "ar1"',model_setup$.__enclos_env__$private$Components))
TempBias <- 'PODat_biasField(main = geometry, model = PODat_bias_field, replicate = year)' %in% model_setup$.__enclos_env__$private$Components

Mesh <- model_setup$.__enclos_env__$private$INLAmesh

rho0 <- Mesh$loc[c(Mesh$segm$bnd$idx[,1], Mesh$segm$bnd$idx[1,1]),1:2] %>% dist(x = .) %>% as.vector(.) %>% max(.)/5
sigma0 <- 1
alpha <- 2; d <- 2; nu <- alpha - d/2
kappa0 <- log(8*nu)/2 - log(rho0)
tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0

spde <- inla.spde2.matern(mesh = Mesh, B.tau = cbind(tau0, nu, -1),
                          B.kappa = cbind(kappa0, -1, 0), constr = FALSE)

numTimes <- length(list.files(paste0(File, 'Data')))

if (!dir.exists(paste0(File, 'ModelRes'))) dir.create(paste0(File, 'ModelRes'))
if (!dir.exists(paste0(File, 'Misc'))) dir.create(paste0(File, 'Misc'))

print(File)

for (i in seq(1, numTimes)) {

print(i)
file <- grep(i, list.files(paste0(File, '/Data/')))
file <- list.files(paste0(File, '/Data/'))[file[1]]
model_setup <- readRDS(paste0(File, '/Data/', file))

numIn <- length(unique(unlist(model_setup$.__enclos_env__$private$temporalVars)))

  if (i == 1) {

    model_setup$.__enclos_env__$private$Components <- sort(model_setup$.__enclos_env__$private$Components)

    model <- fitISDM(model_setup, options = list(bru_max_iter = 1,
                                                 num.threads = 4,
                                                 safe = FALSE,
                                                 inla.mode = 'experimental',
                                                 control.inla = list(
                                                   int.strategy = 'ccd',
                                                   cmin = 0,
                                                   compute.initial.values = TRUE,
                                                   control.vb = list(enable = FALSE),
                                                   strategy = 'adaptive',
                                                   diagonal = 0)))

  } else {

    if (!exists('model')) {

         Mode <- readRDS(paste0(File, 'Misc/ModeE', i-1,'.rds'))
         mode.hyper <- Mode$mode
         names(mode.hyper) <- row.names(Mode)
         Fixed <- readRDS(paste0(File, 'ModelRes/modelRes', i-1,'.rds'))$latent_fixed_effects[[1]]
         fixed.mean <- Fixed$mean
         fixed.prec <- Fixed$sd**-2
         names(fixed.mean) <- names(fixed.prec) <- rownames(Fixed)
         model = list(internal.summary.hyperpar = Mode)

    } else {

    mode.hyper <- model$internal.summary.hyperpar$mode
    names(mode.hyper) <- row.names(model$internal.summary.hyperpar$mode)
    fixed.mean <- model$summary.fixed$mean #%>% list(.)
    fixed.prec <- model$summary.fixed$sd**-2 #%>% list(.)
    names(fixed.mean) <- names(fixed.prec) <- rownames(model$summary.fixed)

    }

    for (cov in names(fixed.mean)[!grepl('intercept', names(fixed.mean))]) {
      model_setup$priorsFixed(cov, mean.linear = fixed.mean[cov], prec.linear = fixed.prec[cov])
    }

    if (ar1Mod) {

          AR1INDP <- grepl('yearTrend', row.names(model$internal.summary.hyperpar)) &  grepl('precision', row.names(model$internal.summary.hyperpar))
          AR1INDR <- grepl('yearTrend', row.names(model$internal.summary.hyperpar)) &  grepl('Rho', row.names(model$internal.summary.hyperpar))

          YEARIND <- grepl('yearTrend\\(',model_setup$.__enclos_env__$private$Components)

          MODEVAR <- tail(model$summary.random$yearTrend$mode, 1)

          model_setup$.__enclos_env__$private$Components[YEARIND] <- 'yearTrend(main = year, model = "ar1", constr = FALSE,hyper = list(prec = norm.internal.sp.ar1P, rho = norm.internal.sp.ar1R, mean = list(initial = MODEVAR)))'

          norm.internal.sp.ar1P <- list(prior = "normal", param = c(model$internal.summary.hyperpar[AR1INDP,1],
                                                                   model$internal.summary.hyperpar[AR1INDP,2]**-2), fixed = FALSE)

          norm.internal.sp.ar1R <- list(prior = "normal", param = c(model$internal.summary.hyperpar[AR1INDR,1],
                                                                    model$internal.summary.hyperpar[AR1INDR,2]**-2), fixed = FALSE)

    }

    ##if both PO and PA
    if (all(c('BTO', 'PODat') %in% model_setup$.__enclos_env__$private$dataSource)) {

    if (any(grepl('_intercept', names(fixed.mean)))) {

    for (int in c('BTO', 'PODat')) {
      model_setup$priorsFixed(Effect = 'Intercept', mean.linear = fixed.mean[paste0(int,'_intercept')], prec.linear = fixed.prec[paste0(int, '_intercept')],datasetName = int)


    }
    }

    BTOCOMP <- grepl('PODat_spatial\\(',model_setup$.__enclos_env__$private$Components)
    model_setup$.__enclos_env__$private$Components[BTOCOMP] <- paste0('PODat_spatial(main = geometry, model = PODat_field, group = year, ngroup = ',numIn,', control.group = list(model = "ar1", hyper = list(rho = list(norm_internal.sp.rho))))')

    BTOIND <- grepl('PODat_spatial', row.names(model$internal.summary.hyperpar)) & grepl('Theta', row.names(model$internal.summary.hyperpar))

    model_setup$spatialFields$datasetFields[[1]] <- inla.spde2.matern(mesh = Mesh,
                                                                     B.tau = cbind(tau0, nu, -1),
                                                                     alpha = 2,
                                                                     B.kappa = cbind(kappa0, -1, 0),
                                                                     constr = TRUE,
                                                                     theta.prior.mean = model$internal.summary.hyperpar[BTOIND,1],
                                                                     theta.prior.prec = model$internal.summary.hyperpar[BTOIND,2]**-2)

    RHOIND <- grepl('PODat_spatial', row.names(model$internal.summary.hyperpar)) & grepl('Group', row.names(model$internal.summary.hyperpar))

    norm_internal.sp.rho <- list(prior = "normal", param = c(model$internal.summary.hyperpar[RHOIND,1],
                                                             model$internal.summary.hyperpar[RHOIND,2]**-2), fixed = FALSE)

    if (length(model_setup$spatialFields$biasFields$PODat) > 1) {

    alpha <- 2; d <- 2; nu <- alpha - d/2
    kappa0 <- log(8*nu)/2 - log(rho0)
    tau0 <- 0.5*(lgamma(nu) - lgamma(nu + d/2) - d/2*log(4*pi)) - log(sigma0) - nu*kappa0

    BIASIND <- grepl('PODat_biasField', row.names(model$internal.summary.hyperpar)) & grepl('Theta', row.names(model$internal.summary.hyperpar))

    model_setup$spatialFields$biasFields$PODat <- inla.spde2.matern(mesh = Mesh,
                                                                    B.tau = cbind(tau0, nu, -1),
                                                                    B.kappa = cbind(kappa0, -1, 0),
                                                                    constr = ifelse(TempBias, FALSE, TRUE),
                                                                    theta.prior.mean = model$internal.summary.hyperpar[BIASIND,1],
                                                                    theta.prior.prec = model$internal.summary.hyperpar[BIASIND,2]**-2)

    }

    }
    else {

        model_setup$priorsFixed(Effect = 'Intercept', mean.linear = fixed.mean[paste0('BTO','_intercept')],
                                prec.linear = fixed.prec[paste0('BTO', '_intercept')],datasetName = 'BTO')

          SHAREDIND <- grepl('shared_spatial', row.names(model$internal.summary.hyperpar)) & grepl('Theta', row.names(model$internal.summary.hyperpar))


          model_setup$spatialFields$sharedField$sharedField <- inla.spde2.matern(mesh = Mesh,
                                                                                 B.tau = cbind(tau0, nu, -1),
                                                                                 B.kappa = cbind(kappa0, -1, 0),
                                                                                 constr = TRUE,
                                                                                 theta.prior.mean = model$internal.summary.hyperpar[SHAREDIND, 1],
                                                                                 theta.prior.prec = model$internal.summary.hyperpar[SHAREDIND, 2]**-2)

          RHOIND <- grepl('shared_spatial', row.names(model$internal.summary.hyperpar)) & grepl('Group', row.names(model$internal.summary.hyperpar))


          norm_internal.sp.rho <- list(prior = "normal", param = c(model$internal.summary.hyperpar[RHOIND,1],
                                                                   model$internal.summary.hyperpar[RHOIND,2]**-2), fixed = FALSE)


          SHAREDCOMP<- grepl('shared_spatial\\(',model_setup$.__enclos_env__$private$Components)
          model_setup$.__enclos_env__$private$Components[SHAREDCOMP] <- paste0('shared_spatial(main = geometry, model = shared_field, group = year, ngroup =', numIn,', control.group = list(model = "ar1", hyper = list(rho = list(norm_internal.sp.rho))))')


    }

    rm(model)
    model_setup$.__enclos_env__$private$Components <- sort(model_setup$.__enclos_env__$private$Components)

    model <- fitISDM(model_setup, options = list(bru_max_iter = 1,
                                                 num.threads = 4,
                                                 safe = FALSE,
                                                 control.mode = list(theta = mode.hyper),
                                                 inla.mode = 'experimental',
                                                 control.inla = list(
                                                   int.strategy = 'ccd',
                                                   cmin = 0,
                                                   control.vb = list(enable = FALSE),
                                                   strategy = 'adaptive',
                                                   diagonal = 0)))

  }


  if (i == numTimes) {

        thetaMode <- model$misc$configs$config[[1]]$theta
        intDesign <- model$joint.hyper
        modNames <- rownames(model$summary.fixed)

        if(!'PODat' %in% model_setup$.__enclos_env__$private$dataSource) {

         spdeMain <- model$bru_info$model$effects$shared_spatial$env$shared_field
         spdeBias <- NULL


        }
        else{

              spdeMain <- model$bru_info$model$effects$BTO_spatial$env$PODat_field

              spdeBias <- model$bru_info$model$effects$PODat_biasField$env$PODat_bias_field
        }


        BLC <- model$summary.spde2.blc
        saveRDS(list(thetaMode = thetaMode,
                     intDesign = intDesign,
                     spdeMain = spdeMain,
                     spdeBias = spdeBias,
                     BLC = BLC), paste0(File, 'finalModEsts.rds'))

        if (!dir.exists(paste0(File, 'FinalRes'))) dir.create(paste0(File, 'FinalRes'))

        if (!is.null(model$bru_info$model$effects[['PODat_spatial']]$main$model)) saveRDS(INLA::inla.spde2.result(inla = model,
                                                                                                                     name = 'PODat_spatial',
                                                                                                                     spde = model$bru_info$model$effects[['PODat_spatial']]$main$model),
                                                                                             file = paste0(File, 'FinalRes/PODat_spatial.rds'))
        else saveRDS(INLA::inla.spde2.result(inla = model,
                                             name = 'shared_spatial',
                                             spde = model$bru_info$model$effects[['shared_spatial']]$main$model),
                     file = paste0(File, 'FinalRes/BTOResult.rds'))


        if (!is.null(model$bru_info$model$effects[['PODat_biasField']]$main$model)) saveRDS(INLA::inla.spde2.result(inla = model,
                                                                                                                       name = 'PODat_biasField',
                                                                                                                       spde = model$bru_info$model$effects[['PODat_biasField']]$main$model),
                                                                                               file = paste0(File, 'FinalRes/BiasResult.rds'))

      marginal_trend <- model$internal.marginals.hyperpar$`Log precision for yearTrend`

      saveRDS(marginal_trend, paste0(File, 'FinalRes/trendResult.rds'))


      saveRDS(model$summary.hyper, paste0(File,'FinalRes/hyperparRes.rds'))


  }

  x <- list(model$mode$x)
  latent_fixed_effects <- list(model$summary.fixed)

  saveRDS(list(x = x, latent_fixed_effects = latent_fixed_effects), paste0(File, 'ModelRes/modelRes',i,'.rds'))
  saveRDS(model$summary.random, paste0(File, 'Misc/RandomE',i,'.rds'))
  saveRDS(model$internal.summary.hyperpar, paste0(File, 'Misc/ModeE',i,'.rds'))
  saveRDS(model$summary.hyperpar, paste0(File, 'Misc/HyperE',i,'.rds'))



  print(paste("Model", i, "is done."))
  print(File)


}

