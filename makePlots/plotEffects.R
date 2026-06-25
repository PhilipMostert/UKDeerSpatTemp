library(ggplot2)
Folder <- 'xx/'
yearsIn <-  2005:2024
Mesh <- readRDS('Data/MeshCoarse.rds')

Random <- list()
POTrend <- list()
MatCorrelation <- list()
Range <- list()
Fixed <- list()

spec <- c('Capreoluscapreolus', 'Cervuselaphus', 'Damadama', 'Muntiacusreevesi')

whichSpat <- 'PODat_spatial'#'BTOResult'
distRange <- seq(0, 100)

for (sp in spec) {

  nIn <- length(list.files(paste0(Folder, spec, '/modelRes/')))

  specPlot <- switch(sp,
                     Muntiacusreevesi = 'Muntiacus reevesi',
                     Damadama = 'Dama dama',
                     Capreoluscapreolus = 'Capreolus capreolus',
                     Cervuselaphus = 'Cervus elaphus')

  R <- list()
  for (l in 1:nIn) {

    R[[l]] <- readRDS(paste0(Folder, sp, '/Summaries/sumRandom', l, '.rds'))$yearTrend

    if (l == nIn) int <- readRDS(paste0(Folder, sp, '/Summaries/fixedRes', l, '.rds'))['BTO_intercept',]

  }

  R <- do.call(rbind, R)

  R$yearsIn <- yearsIn
  R$Species <- specPlot

  Random[[sp]] <- R

  MC <- readRDS(paste0(Folder, sp, '/Summaries/', whichSpat, '.rds'))

  out <- inlabru::materncov.bands(manifold = Mesh$Mesh,
                                  dist = distRange,
                                  log.range = list(mean = MC$summary.log.range.nominal$mean,
                                                   sd = MC$summary.log.range.nominal$sd))

  MatCorrelation[[sp]] <- data.frame(x = distRange, q0.5 = out$median, lower = out$lower, upper = out$upper, median = out$median, Species = specPlot)

  RangeMarg <- MC$marginals.range.nominal[[1]]
  med <- INLA::inla.qmarginal(0.5, RangeMarg)

  inner.x <- seq(0, 80, length.out = 1000)
  Range[[sp]] <- data.frame(
    x = inner.x,
    y = INLA::inla.dmarginal(inner.x, RangeMarg), Species = specPlot
  )

  FixedE <- readRDS(paste0(Folder, sp, '/Summaries/fixedRes', nIn, '.rds'))
  FixedE <- FixedE[!grepl('intercept',row.names(FixedE)),]

  FixedE$Species <- specPlot
  FixedE$Variable <- row.names(FixedE)

  FixedE$Variable[FixedE$Variable == 'SGrassland'] <- 'Semi-natural grassland'
  FixedE$Variable[FixedE$Variable == 'time_total'] <- 'Duration'
  FixedE$Variable[FixedE$Variable == 'transects_year'] <- 'Number of visits'

  FixedE$Type <- 'Observation'
  FixedE$Type[FixedE$Variable != 'Duration'] <- 'Ecological'
  Fixed[[sp]] <- FixedE[order(row.names(FixedE)),]

}

Random <- do.call(rbind, Random)
MatCorrelation <- do.call(rbind, MatCorrelation)
Range <- do.call(rbind, Range)
Fixed <- do.call(rbind, Fixed)

source('Misc/theme.R')

yearTheme <- themePlot
#yearTheme$legend.position <- 'none'
yearTheme$axis.text.x$angle <- 45
yearTheme$axis.text.x$hjust <- 1
#yearTheme$axis.title.y <- element_text(vjust  = 1)
yearTheme$legend.position <- 'bottom'
yearTheme$legend.justification <- 'center'
yearTheme$legend.text <- element_text(size = 12*4)
yearTheme$legend.title <- element_text(size = 12*5)
yearTheme$strip.background <- element_blank()

##plot year trend
yearPlot <- ggplot(data = Random, aes(x = yearsIn, y = `0.5quant`)) +
  geom_point(size = 0.75, aes(col = Species)) +
  scale_x_continuous(limits = c(min(yearsIn),max(yearsIn)+0.25), expand = c(0, 0),
                     labels = as.character(yearsIn)[seq(1, length(yearsIn), by = 3)],
                     breaks = yearsIn[seq(1, length(yearsIn), by = 3)]) +
  scale_colour_manual(values = c('#5FAD56', '#F78154', '#B4436C', '#067BC2')) +
  ylab(expression(tilde(gamma[t]))) +
  xlab('Year') +
  labs(tag = 'B)') +
  geom_line(aes(col = Species), size = 1.2) +
  geom_ribbon(aes(ymin=  `0.025quant`, ymax= `0.975quant`, fill = Species), alpha = 0.2) +
  yearTheme +
  geom_vline(xintercept = c(2005, 2009, 2014, 2019, 2024), linetype = 'dashed', alpha = 0.5) +
  #geom_vline(xintercept = c(2005, 2008, 2012, 2016, 2020, 2024), linetype = 'dashed', alpha = 0.5) +
  #geom_vline(xintercept = c(2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023), linetype = 'dashed', alpha = 0.5) +
  geom_hline(yintercept  = 0, linetype = 'dashed') +
  facet_wrap(~Species) +
  guides(colour = guide_legend(override.aes = list(linetype = "blank", size = 10, shape = 15)))

#Mat cor
marCorPlot <- ggplot() +
  geom_line(data = MatCorrelation, aes(x = x, y = median, col = Species)) +
  geom_ribbon(data = MatCorrelation, aes(x = x, y = median, ymin= lower, ymax= upper, fill = Species), alpha = 0.1) +
  scale_x_continuous(limits = c(0, max(distRange)), expand = c(0, 0)) +
  scale_colour_manual(values = c('#5FAD56', '#F78154', '#B4436C', '#067BC2')) +
  scale_fill_manual(values = c('#5FAD56', '#F78154', '#B4436C', '#067BC2')) +
  xlab('Distance (km)') + ylab('Correlation') +
  themePlot + ggtitle(expression(Matérn~correlation~`for`~zeta(s,t))) +
  guides(fill = guide_legend(override.aes = list(linetype = "blank", size = 10, shape = 1, alpha = 1)))

##Range
rangePlot <- ggplot() +
  geom_line(data = Range, aes(x = x, y = y, col = Species)) +
  geom_ribbon(data = Range, aes(x = x, ymin = 0, ymax = y, fill = Species), alpha = 0.1) +
  scale_x_continuous(limits = c(8, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.27),expand = c(0,0)) +
  scale_x_continuous(limits = c(6, 80 + 10), expand = c(0, 0)) +
  scale_colour_manual(values = c('#5FAD56', '#F78154', '#B4436C', '#067BC2')) +
  scale_fill_manual(values = c('#5FAD56', '#F78154', '#B4436C', '#067BC2')) +
  xlab('Distance (km)') + ylab('Density') + ggtitle(expression(Range~estimates~of~zeta(s,t))) +
  themePlot + guides(fill = guide_legend(override.aes = list(linetype = "blank", size = 10, shape = 1, alpha = 1)))

Fixed$Variable <- factor(Fixed$Variable, levels = c(unique(Fixed$Variable[Fixed$Variable != 'Duration']), 'Duration'), ordered = TRUE)
#levels(Fixed$Variable) <- c(levels(Fixed$Variable)[levels(Fixed$Variable) != 'Duration'], 'Duration')
#Fixed

themePlot$legend.position <- 'none'
fixedPlot <- ggplot() +
  geom_point(data = Fixed, aes(x = `0.5quant`, y = Variable, col = Species), position=position_dodge(width=0.2), size = 4) +
  geom_errorbar(data = Fixed, aes(x = `0.5quant`, y = Variable, xmin = `0.025quant`, xmax = `0.975quant`, col = Species),
  position=position_dodge(width=0.2), width = 0, size = 1.1) +
  scale_colour_manual(values = c('#5FAD56', '#F78154', '#B4436C', '#067BC2'), ) + ##067BC2
  geom_vline(xintercept = 0, lty = 2) +
  scale_y_discrete(limits = sort(Fixed$Variable)) +
  ylab('Fixed effect') + xlab('Median') +
  geom_vline(xintercept = 0, lty = 2) +
  #xlim(-1.5, 1.5) +
  themePlot +
   labs(tag = 'A)') #+
  # guides(colour = guide_legend(override.aes = list(linetype = "blank", size = 5, shape = 15)))


ggsave(filename = paste0(Folder, 'Plots/FixedEffects.png'), dpi = 300, plot = fixedPlot, scale = 1.2, width = 3618, height = 2492, units = 'px')
ggsave(filename = paste0(Folder, 'Plots/yearPlot.png'), dpi = 300, plot = yearPlot, scale = 1.2, width = 3618, height = 2492, units = 'px')
