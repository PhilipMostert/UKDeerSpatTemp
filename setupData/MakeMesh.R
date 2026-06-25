library(sf)
library(terra)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(ggplot2)
library(fmesher)
library(inlabru)
crs <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
UK <- giscoR::gisco_get_countries(country = 'United Kingdom', resolution = '10')
UK <- rmapshaper::ms_filter_islands(UK, min_area = 100000000000)

UKsimp <-  rmapshaper::ms_simplify(UK, keep = 0.5)
UKsimp <- st_transform(UKsimp, crs)
UKsimp <- st_buffer(UKsimp, 5)

max.edge = 16

Meshpars <- list(max.edge = c(max.edge*1.6, max.edge*4), #*4 #0.75 #1.5
                 offset = c(max.edge, max.edge*10), #85
                 cutoff = max.edge/2) #2

Mesh <- fmesher:::fm_mesh_2d(boundary = UKsimp,
                             max.edge = Meshpars$max.edge,max.n.strict = 1500,
                             offset = Meshpars$offset, cutoff = Meshpars$cutoff, crs = crs, min.angle = 20)

IPS <- fm_int(fm_subdivide(Mesh, 8), samplers = UKsimp)
saveRDS(list(UK = UKsimp, Mesh = Mesh, IPS = IPS), 'Data/MeshCoarse.rds')
