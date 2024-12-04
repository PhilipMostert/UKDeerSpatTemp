library(sf)
library(terra)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(ggplot2)
library(inlabru)
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)
UKsimp <- st_buffer(UKsimp, 5)
# max.edge = 30#20 #12
# Meshpars <- list(max.edge = c(max.edge*0.75, max.edge*1.25), #*4
#                  offset = c(max.edge, max.edge*2.5), #85
#                  cutoff = max.edge/2) #5

max.edge = 16#20 #12
Meshpars <- list(max.edge = c(max.edge*0.9, max.edge*4), #*4 #0.75 #1.5
                 offset = c(max.edge, max.edge*8), #85
                 cutoff = max.edge/3) #3

Mesh <- fmesher:::fm_mesh_2d(boundary = UKsimp, max.edge = Meshpars$max.edge,
                            offset = Meshpars$offset, cutoff = Meshpars$cutoff, crs = crs)

saveRDS(list(UK = UKsimp, Mesh = Mesh), 'Data/MeshFine.rds')

Meshpars <- list(max.edge = c(max.edge*2, max.edge*5), #*4 #0.75 #1.5
                 offset = c(max.edge, max.edge*10), #85
                 cutoff = max.edge/1.8) #2

Mesh <- fmesher:::fm_mesh_2d(boundary = UKsimp, max.edge = Meshpars$max.edge,max.n.strict = 1500,
                             offset = Meshpars$offset, cutoff = Meshpars$cutoff, crs = crs)

saveRDS(list(UK = UKsimp, Mesh = Mesh), 'Data/MeshCoarse.rds')

Meshpars <- list(max.edge = c(max.edge*0.5, max.edge*2), #*4 #0.75 #1.5
                 offset = c(max.edge, max.edge*5), #85
                 cutoff = max.edge/3) #3

Mesh <- fmesher:::fm_mesh_2d(boundary = UKsimp, max.edge = Meshpars$max.edge, min.angle = 26,
                             offset = Meshpars$offset, cutoff = Meshpars$cutoff, crs = crs)

saveRDS(list(UK = UKsimp, Mesh = Mesh), 'Data/MeshUFine.rds')
# ggplot() + gg(fm_transform(Mesh, 'WGS84')) +
#   coord_sf(crs = 'WGS84') + theme_bw() +
#   ggtitle('Triangulated mesh used to \nestimate the spatial effect') +
#   xlab('Lon') + ylab('Lat') +
#   theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text = element_text(size = 10), axis.title = element_text(size = 15))
