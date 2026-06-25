library(showtext)
library(ggplot2)

Folder <- 'xx/'

Files <- c('Full', 'IID', 'noBias', 'IIDnoBias', 'BTO', 'BTOIID')

species <- c('Capreoluscapreolus', 'Cervuselaphus', 'Damadama', 'Muntiacusreevesi')

Results <- list()

for (fl in Files) {

  modUsed <- switch(fl,
                    Full = 'AR1 with bias field',
                    IID = 'IID with bias field',
                    noBias = 'AR1',
                    IIDnoBias = 'IID',
                    BTO = 'BTO AR1',
                    BTOIID = 'BTO IID')

  for (spec in gsub(' ', '', species)) {

        specUsed <- switch(spec,
                           Capreoluscapreolus = 'Capreolus capreolus',
                           Cervuselaphus = 'Cervus elaphus',
                           Damadama = 'Dama dama',
                           Muntiacusreevesi = 'Muntiacus reevesi')

    if (file.exists(paste0(Folder, fl, '/', spec, '/intensitySamps.rds'))) {


    res <- readRDS(paste0(Folder, fl, '/', spec, '/intensitySamps.rds'))

    year <- as.factor(res$year)

    pres <- res$pres

    dev <- apply(res[, grepl('Sample:', colnames(res))], MARGIN = 2, FUN = function(x) {


      deviance(glm(pres ~ x , family = poisson))

    })

    } else dev <- data.frame(meanD = NA)

    Results[[fl]][[spec]] <-  data.frame(meanD = dev, Species = specUsed, Model = modUsed)


  }

}

plotData <- do.call(rbind, unlist(Results, recursive = FALSE))
plotData$Model <- forcats::fct_reorder(plotData$Model, nchar(plotData$Model), .desc = TRUE)


showtext_auto()

f1 <- "Open sans"
font_add_google(f1, f1)

ggplot(data = plotData, aes(x = meanD, y = Model, fill = Model)) +
  geom_boxplot(outlier.size = 0.001) +
  scale_fill_manual(values = c('#D5BBB1', '#6E9075', '#80A4ED', '#E07A5F', '#F2CC8F', '#EF767A')) +
  xlab('Deviance') +
  facet_wrap(~Species, scales = 'free') +
  theme_minimal() +
  theme(text = element_text(family = f1),
        strip.text= element_text(size = 15*4),
        axis.text.x = element_text(size = 12*3, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12*3),
        axis.title = element_text(size = 15*4),
        legend.position = "none",
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.justification='center')

ggsave(paste0(Folder, '/CrossValidation.png'), scale = 1)


