library(ggplot2)
library(terra)
library(sf)
library(inlabru)
library(tidyterra)
library(showtext)
showtext_auto()
f1 <- "Open sans"
font_add_google(f1, f1)

Folder <- 'xx/'
spec <- c('Capreoluscapreolus', 'Cervuselaphus', 'Damadama', 'Muntiacusreevesi')

UK <- st_read('Data/NUTS1_Jan_2018_UGCB_in_the_UK_2022_-7602393063213847322/')
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"

UK <- st_transform(UK, crs)
UK <- UK[UK$nuts118nm != 'Northern Ireland',]
UK <- rmapshaper::ms_filter_islands(UK, min_area = 1000)

UKsimp <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UKsimp <- rmapshaper::ms_filter_islands(UKsimp, min_area = 100000000000)
UKsimp <- st_transform(UKsimp, crs = crs)

UK <- rmapshaper::ms_clip(UK, UKsimp)

intensityMedian <- list()
intensitySD <- list()

biasMedian <- list()
biasSD <- list()

diffMedian <- terra::rast()
diffSD <- terra::rast()

intensity2021 <- terra::rast()
sd2021 <- terra::rast()

themeChange <- theme_void() +
  theme(text = element_text(family = f1),
        panel.background=element_blank(),
        legend.key.size = unit(0.85, 'cm'),
        strip.text= element_text(size = 15*1.7),
        axis.text.x = element_text(size = 12*1.5),
        axis.text.y = element_text(size = 12*1.5),
        axis.title = element_text(size = 15*1.5),
        legend.text = element_text(size = 12*1.8),
        legend.title = element_text(size = 15*2))

themeInt <- themeChange
themeInt$plot.title <-  element_text(hjust = 0.5, size = 15*4, vjust = 1) # vjust = 2
themeInt$strip.background <- element_blank()
themeInt$legend.title <- element_text(size = 15*3)
themeInt$axis.text.x <- themeInt$axis.text.y <-  element_text(size = 12*2.2) #12*1.5
themeInt$strip.text <- element_text(size = 15*2.5) #15*2
themeInt$legend.justification <- c(0.0, 0.5)
themeInt$legend.margin <-  margin(t = 0, unit = "pt")
plot_ratio <- tmaptools::get_asp_ratio(UK)

themeSD <- themeInt
themeSD$plot.title <- element_blank()

for (sp in spec) {

  specPlot <- switch(sp,
                     Muntiacusreevesi = 'Muntiacus reevesi',
                     Damadama = 'Dama dama',
                     Capreoluscapreolus = 'Capreolus capreolus',
                     Cervuselaphus = 'Cervus elaphus')

  intensityMedian[[specPlot]] <- readRDS(paste0(Folder, sp, '/Results/', sp, '_linPredMedian.rds'))
  intensitySD[[specPlot]] <- readRDS(paste0(Folder, sp, '/Results/', sp, '_linPredSD.rds'))

  if (file.exists(paste0(Folder, sp, '/Results/', sp ,'_biasMedian.rds'))) {

  biasMedian[[specPlot]] <- readRDS(paste0(Folder, sp, '/Results/', sp ,'_biasMedian.rds'))
  biasSD[[specPlot]] <- readRDS(paste0(Folder, sp, '/Results/', sp ,'_biasSD.rds'))

  }


  medDiff <- subset(intensityMedian[[specPlot]], nlyr(intensityMedian[[specPlot]])) - subset(intensityMedian[[specPlot]], 1)
  sdDiff <- subset(intensitySD[[specPlot]], nlyr(intensitySD[[specPlot]])) - subset(intensitySD[[specPlot]], 1)

  names(medDiff) <- names(sdDiff) <- specPlot
  diffMedian <- c(diffMedian, medDiff)
  diffSD <- c(diffSD, sdDiff)

  int2021 <- subset(intensityMedian[[specPlot]], nlyr(intensityMedian[[specPlot]]))
  s2021 <- subset(intensitySD[[specPlot]], nlyr(intensitySD[[specPlot]]))

  names(int2021) <- names(s2021) <- specPlot

  intensity2021 <- c(intensity2021, int2021)
  sd2021 <- c(sd2021, s2021)

  intMedPlot <- ggplot() +
    geom_spatraster(data = intensityMedian[[specPlot]]) +
    facet_wrap(~lyr, nrow = 2) +
    gg(st_boundary(UK), lwd = 1) +
    labs(fill = 'Median') +
    viridis::scale_fill_viridis(option="magma", na.value = NA) +
    ggtitle(paste0('Log intensity for ', specPlot)) +
    themeInt

  intSDPlot <- ggplot() +
    geom_spatraster(data = intensitySD[[specPlot]]) +
    facet_wrap(~lyr, nrow = 2) +
    gg(st_boundary(UK), lwd = 1) +
    labs(fill = 'SD') +
    viridis::scale_fill_viridis(option="magma", na.value = NA) +
    ggtitle(paste0('Log intensity for ', specPlot)) +
    themeSD

  cowplot::set_null_device("png")
  cowplot::plot_grid(intMedPlot, NULL ,intSDPlot, rel_heights = c(1, -0.1, 1), ncol = 1, align = 'hv')

  ggsave(filename = paste0(Folder, 'Plots/',sp,'/', sp, '_intensity.png'),
         height = plot_ratio*50, width = 50/1.4, units = 'cm')

  if (length(biasMedian) > 0) {

  biasMedPlot <- ggplot() +
    geom_spatraster(data = biasMedian[[specPlot]]) +
    facet_wrap(~lyr, nrow = 2) +
    gg(st_boundary(UK), lwd = 1) +
    labs(fill = 'Median') +
    viridis::scale_fill_viridis(option="magma", na.value = NA) +
    ggtitle(paste0('Log intensity for ', specPlot)) +
    themeInt

  biasSDPlot <- ggplot() +
    geom_spatraster(data = biasSD[[specPlot]]) +
    facet_wrap(~lyr, nrow = 2) +
    gg(st_boundary(UK), lwd = 1) +
    labs(fill = 'SD') +
    viridis::scale_fill_viridis(option="magma", na.value = NA) +
    ggtitle(paste0('Log intensity for ', specPlot)) +
    themeSD
    #themeInt

  cowplot::set_null_device("png")
  cowplot::plot_grid(biasMedPlot, NULL ,biasSDPlot, rel_heights = c(1, -0.1, 1), ncol = 1, align = 'hv')

  ggsave(filename = paste0(Folder, 'Plots/',sp,'/', sp, '_bias.png'),
         height = plot_ratio*50, width = 50/1.4, units = 'cm')

  }

}

#Blue vs orange is colourblind friendly
#Light Blue & Pink

coloursDiff = c('#71092C','#BD0F49','#EE2F6F','#FAC7D8', '#FAC7D8','#FFFFFF','#C2FFF0' ,'#70FFDB', '#1FFFC7', '#00CC99', '#007A5C')
DiffPlot <- ggplot() +
  geom_spatraster(data = diffMedian) +
  facet_wrap(~lyr, nrow = 1) +
  scale_fill_steps2(low = coloursDiff[1], high = coloursDiff[11], midpoint = 0, mid = '#FFFFFF', n.breaks = 10, na.value = "transparent") +
  gg(st_boundary(UK), lwd = 0.25) +
  gg(st_boundary(UKsimp)) +
  labs(fill = expression(Delta~"Median")) +
  themeChange


DiffSD <- ggplot() +
  geom_spatraster(data = diffSD) +
  facet_wrap(~lyr, nrow = 1) +
  scale_fill_steps2(low = coloursDiff[1], high = coloursDiff[11], midpoint = 0, mid = '#FFFFFF', n.breaks = 10, na.value = "transparent") +
  gg(st_boundary(UK), lwd = 0.25) +
  gg(st_boundary(UKsimp)) +
  labs(fill = expression(Delta~"SD"~phantom(xxx))) +
  themeChange

cowplot::set_null_device("png")
cowplot::plot_grid(DiffPlot, DiffSD, ncol = 1, align = 'hv')

ggsave(paste0(Folder,'/Plots/Changes.png'), dpi = 300, height = 2077, width = 3015, units = 'px')

