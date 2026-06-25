Folder <- 'xx/'

speciesIn <- c('Capreoluscapreolus', 'Cervuselaphus', 'Damadama', 'Muntiacusreevesi')

Mesh <- readRDS('Data/MeshCoarse.rds')

Results <- list()

whichMax <- list.files(paste0(Folder, speciesIn[1], '/Data'))
whichMax <- max(as.numeric(gsub("[^\\d]+", "", whichMax, perl=TRUE)),na.rm = TRUE)

for (species in speciesIn) {

  print(species)

  Hyper <- readRDS(paste0(Folder, species, '/FinalRes/hyperparRes.rds'))

  spatSum <- readRDS(paste0(Folder, species, '/FinalRes/PODat_spatial.rds'))

  print('Results Spatial:')
  spatRange <- INLA:::inla.zmarginal(marginal = spatSum$marginals.range.nominal[[1]])
  spatSD <- INLA:::inla.zmarginal(marginal = spatSum$marginals.variance.nominal[[1]])

  if (file.exists(paste0(Folder, species, '/FinalRes/BiasResult.rds'))) {

    print('Results Bias:')
    spatSum <- readRDS(paste0(Folder, species, '/FinalRes/BiasResult.rds'))
    BiasRange <- INLA:::inla.zmarginal(marginal = spatSum$marginals.range.nominal[[1]])
    BiasSD <- INLA:::inla.zmarginal(marginal = spatSum$marginals.variance.nominal[[1]])

  }

  print('Hyperparameters:')
  print(Hyper[!grepl('Theta', row.names(Hyper)),])

  print('Variance trend:')
  precTrend <- readRDS(paste0(Folder, species, '/FinalRes/trendResult.rds'))
  varTrend <- INLA:::inla.tmarginal(function(x) 1/exp(x), precTrend)
  INLA:::inla.zmarginal(varTrend)

  cat('\n\n')


}

