args <- commandArgs(trailingOnly = TRUE)

Folder <- args[1]
print(Folder)

FolderUsed <- paste0(Folder, 'ModelRes/')

modelEsts <- list(x = list(), latent_fixed_effects = list())

numEsts <- length(grepl('ModelRes', list.files((FolderUsed))))

for (i in 1:numEsts) {

  ests <- readRDS(paste0(FolderUsed, 'modelRes', i, '.rds'))
  modelEsts$x[[i]] <- ests$x[[1]]
  modelEsts$latent_fixed_effects[[i]] <- ests$latent_fixed_effects[[1]]

}

saveRDS(modelEsts, paste0(Folder, 'modEsts.rds'))


rm(list = ls())
