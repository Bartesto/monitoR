#### Code used to create and generate internal data for both rivers
library(raster)
library(ggplot2)
library(gstat)
library(rgdal)
library(rgeos)
library(scales)
library(tidyverse)
library(janitor)
#library(formatR)
library(grid)
library(gridExtra)
library(readxl)
library(metR)
library(ggpubr)
library(lubridate)
library(tidyxl)
library(unpivotr)
library(splancs)

#####Swan River
## Create site data - locations, depths and distances
sites <- readOGR(dsn = "../vectors/swan100",
                 layer = "SwanDepthDistProfile_MGA50_100m",
                 stringsAsFactors = FALSE)
S_sitesdf <- as.data.frame(sites@data) %>%
  dplyr::mutate(site = ifelse(site == "BWR10", "BWR", site)) %>%
  dplyr::arrange(dist_mouth)

## Create oxygen location sites
S_oxy_locs <- S_sitesdf %>%
  dplyr::filter(site == "VIT" | site == "CAV") %>%
  dplyr::mutate(x = dist_mouth/1000,
                y = -6) %>%
  dplyr::select(x, y)

## Bottoms for plotting
S_bottom <- data.frame(x = c(-1, S_sitesdf$dist_mouth/1000,
                             51.6, 51.6, -1),
                     y = c(-10 ,S_sitesdf$adj_depth,
                           -1.9, -22.1, -22.1),
                     id = 1)

S_bottom_nar <- data.frame(x = c(20.95 , S_sitesdf[S_sitesdf$dist_mouth/1000 >= 21, 8]/1000,
                                 51.6, 51.6, 20.95),
                         y = c(-3.9 , S_sitesdf[S_sitesdf$dist_mouth/1000 >= 21, 9],
                               -1.9, -10.05, -10.05),
                         id = 1)
## Irregular grid creation
# for all river
b_all <- data.frame(x = c(-1.1, S_sitesdf$dist_mouth/1000,
                          51.6, 51.6, -1.1),
                    y = c(-10 ,S_sitesdf$adj_depth,
                          -1.9, -22.2, -22.2)) #base and sides adjusted so buffer doesn't impact

# convert to spatial poly to allow buffering
p_all = Polygon(b_all)
ps_all = Polygons(list(p_all),1)
sps_all = SpatialPolygons(list(ps_all))

# buffering in rgeos to allow for negative buffer
buff_sps_all <- rgeos::gBuffer(sps_all, width = -0.1) # neg buffer

# extract coords and convert to df
buff_spdf_all <- as.data.frame(buff_sps_all@polygons[[1]]@Polygons[[1]]@coords)

# for above Narrows
b_nar <- data.frame(x = c(20.95 , S_sitesdf[S_sitesdf$dist_mouth/1000 >= 21, 8]/1000,
                          51.6, 51.6, 20.95),
                    y = c(-3.9 , S_sitesdf[S_sitesdf$dist_mouth/1000 >= 21, 9],
                          -1.9, -10.05, -10.05))

# convert to spatial poly to allow buffering
p_nar = Polygon(b_nar)
ps_nar = Polygons(list(p_nar),1)
sps_nar = SpatialPolygons(list(ps_nar))

# buffering in rgeos to allow for negative buffer
buff_sps_nar <- rgeos::gBuffer(sps_nar, width = -0.1) # neg buffer

# extract coords and convert to df
buff_spdf_nar <- as.data.frame(buff_sps_nar@polygons[[1]]@Polygons[[1]]@coords)

## Create template grids for interpolations
x_range <- as.numeric(c(-1, 51.5))
y_range <- as.numeric(c(-22, 0))
x_range_nar <- as.numeric(c(21, 51.5))
y_range_nar <- as.numeric(c(-10, 0))


# create an empty regular grid of values for whole of river
grd1 <- expand.grid(x = seq(from = x_range[1],
                            to = x_range[2],
                            by = 0.1),
                    y = seq(from = y_range[1],
                            to = y_range[2],
                            by = 0.1))  # expand points to grid

# create an empty regular grid of values for narrows and up
grd2 <- expand.grid(x = seq(from = x_range_nar[1],
                            to = x_range_nar[2],
                            by = 0.1),
                    y = seq(from = y_range_nar[1],
                            to = y_range_nar[2],
                            by = 0.1))  # expand points to grid

# create coord lists
b_list_all <- list(x = buff_spdf_all$x, y = buff_spdf_all$y)
b_list_nar <- list(x = buff_spdf_nar$x, y = buff_spdf_nar$y)

# create irregular grid following buffered bottom - all river
S_grd_all <- grd1[!splancs::inout(grd1, b_list_all),]
coordinates(S_grd_all) <- ~x + y
gridded(S_grd_all) <- TRUE

# create irregular grid following buffered bottom - from Narrows
S_grd_nar <- grd2[!splancs::inout(grd2, b_list_nar),]
coordinates(S_grd_nar) <- ~x + y
gridded(S_grd_nar) <- TRUE

## Reclassification matrices for binning interpolated rasters
# salinity
aSal <- seq(0, 42, 2)
bSal <- rep(aSal, each = 3)
reclass_dfSal <- c(-1, bSal, 44, 44)
reclass_mSal <- matrix(reclass_dfSal,
                       ncol = 3,
                       byrow = TRUE)
# dissolved oxygen
aDo <- seq(0, 17, 1)
bDo <- rep(aDo, each = 3)
reclass_dfDo <- c(-1, bDo, 18, 18)
reclass_mDo <- matrix(reclass_dfDo,
                      ncol = 3,
                      byrow = TRUE)
# temperature
aT <- seq(0, 33, 1)
bT <- rep(aT, each = 3)
reclass_dfT <- c(-1, bT, 34, 34)
reclass_mT <- matrix(reclass_dfT,
                     ncol = 3,
                     byrow = TRUE)
# chlorophyll
aC <- seq(20, 80, 20)
bC <- rep(aC, each = 3)
reclass_dfC <- c(0, bC, 120, 120, 120, 200, 200, 200, 1000, 1000)
reclass_mChl <- matrix(reclass_dfC,
                       ncol = 3,
                       byrow = TRUE)

reclass_matrices <- list(reclass_mSal = reclass_mSal,
                         reclass_mDo = reclass_mDo,
                         reclass_mT = reclass_mT,
                         reclass_mChl = reclass_mChl)

## Create colour breaks for metrics
sal_brk <- as.character(seq(2, 42, 2))
do_mg_l_brk <- as.character(seq(1, 17, 1))
chl_brk <- c(as.character(seq(20, 80, 20)), "120", "200", "1000")
temp_brk <- as.character(seq(11, 33, 1))


#####Canning River
## Create site data - locations, depths and distances
csites <- readOGR(dsn = "../vectors/canning100",
                  layer = "CanningDepthDistProfile_MGA50_100m",
                  stringsAsFactors = FALSE)
C_sitesdf <- as.data.frame(csites@data) %>%
  dplyr::arrange(dist_bridg)

## Create oxygen location sites
C_oxy_locs <- C_sitesdf %>%
  dplyr::filter(site == "BAC" | site == "NIC") %>%
  dplyr::mutate(x = dist_bridg/1000,
                y = -6.8) %>%
  dplyr::select(x, y)

## Bottoms for plotting
C_bottom_weir <- data.frame(x = c(C_sitesdf$dist_bridg[C_sitesdf$dist_bridg >= 500 & C_sitesdf$dist_bridg <= 11100]/1000,
                                  11.2, 11.2, 11.3, 11.3,
                                  C_sitesdf$dist_bridg[C_sitesdf$dist_bridg > 11300]/1000,
                                  15.95, 15.95, 0.5),
                            y = c(C_sitesdf$adj_depth[C_sitesdf$dist_bridg >= 500 & C_sitesdf$dist_bridg <= 11100],
                                  -1.81, 0.1, 0.1, -1.51,
                                  C_sitesdf$adj_depth[C_sitesdf$dist_bridg > 11300],
                                  -1, -7.1, -7.1),
                            id = 1)

C_bottom_open <- data.frame(x = c(C_sitesdf$dist_bridg[C_sitesdf$dist_bridg >= 500 & C_sitesdf$dist_bridg <= 11100]/1000,
                                  11.2, 11.2, 11.3, 11.3,
                                  C_sitesdf$dist_bridg[C_sitesdf$dist_bridg >= 11400]/1000,
                                  15.95, 15.95, 0.5),
                            y = c(C_sitesdf$adj_depth[C_sitesdf$dist_bridg >= 500 & C_sitesdf$dist_bridg <= 11100],
                                  -1.81, -0.5, -0.5, -1.51,
                                  C_sitesdf$adj_depth[C_sitesdf$dist_bridg >= 11400],
                                  -1, -7.1, -7.1),
                            id = 1)

## Irregular grid creation
# below weir
b_low <- data.frame(x = c(0.4, C_sitesdf$dist_bridg[C_sitesdf$dist_bridg >= 500 & C_sitesdf$dist_bridg <= 11300]/1000,
                          11.3, 0.4),
                    y = c(-4.76, C_sitesdf$adj_depth[C_sitesdf$dist_bridg >= 500 & C_sitesdf$dist_bridg <= 11300],
                          -7.1, -7.1)) #base and sides adjusted so buffer doesn't impact

# convert to spatial poly to allow buffering
p_low = Polygon(b_low)
ps_low = Polygons(list(p_low),1)
sps_low = SpatialPolygons(list(ps_low))

# buffering in rgeos to allow for negative buffer
buff_sps_low <- rgeos::gBuffer(sps_low, width = -0.2) # neg buffer

# extract coords and convert to df
buff_spdf_low <- as.data.frame(buff_sps_low@polygons[[1]]@Polygons[[1]]@coords)

# above weir
b_upper <- data.frame(x = c(C_sitesdf$dist_bridg[C_sitesdf$dist_bridg >= 11200]/1000,
                            15.70524, 11.2),
                      y = c(C_sitesdf$adj_depth[C_sitesdf$dist_bridg >= 11200],
                            -7.1, -7.1)) #base and sides adjusted so buffer doesn't impact

# convert to spatial poly to allow buffering
p_up = Polygon(b_upper)
ps_up = Polygons(list(p_up),1)
sps_up = SpatialPolygons(list(ps_up))

# buffering in rgeos to allow for negative buffer
buff_sps_up <- rgeos::gBuffer(sps_up, width = -0.2) # neg buffer

# extract coords and convert to df
buff_spdf_up <- as.data.frame(buff_sps_up@polygons[[1]]@Polygons[[1]]@coords)

# whole river
b_all <- data.frame(x = c(0.4, C_sitesdf$dist_bridg[C_sitesdf$dist_bridg >= 500]/1000,
                          15.95, 15.95, 0.4),
                    y = c(-4.76, C_sitesdf$adj_depth[C_sitesdf$dist_bridg >= 500],
                          -1, -7.1, -7.1)) #base and sides adjusted so buffer doesn't impact

# convert to spatial poly to allow buffering
p_all = Polygon(b_all)
ps_all = Polygons(list(p_all),1)
sps_all = SpatialPolygons(list(ps_all))

# buffering in rgeos to allow for negative buffer
buff_sps_all <- rgeos::gBuffer(sps_all, width = -0.2) # neg buffer

# extract coords and convert to df
buff_spdf_all <- as.data.frame(buff_sps_all@polygons[[1]]@Polygons[[1]]@coords)

# 1# lower canning
# establish an extent within which you want to interpolate
lc_x_range <- as.numeric(c(0.5, 11.2))  # min/max x of the interpolation area
lc_y_range <- as.numeric(c(-7, 0))  # min/max y of the interpolation area

# 2# upper canning
# establish an extent within which you want to interpolate
uc_x_range <- as.numeric(c(11.3, 15.95))  # min/max x of the interpolation area
uc_y_range <- as.numeric(c(-7, 0))  # min/max y of the interpolation area

# 3# whole of river for open weir
# establish an extent within which you want to interpolate
all_x_range <- as.numeric(c(0.5, 15.95))  # min/max x of the interpolation area
all_y_range <- as.numeric(c(-7, 0))  #


# create an empty grid for lower canning
grd1 <- expand.grid(x = seq(from = lc_x_range[1],
                            to = lc_x_range[2],
                            by = 0.02),
                    y = seq(from = lc_y_range[1],
                            to = lc_y_range[2],
                            by = 0.1))  # expand points to grid

# create an empty grid for upper canning
grd2 <- expand.grid(x = seq(from = uc_x_range[1],
                            to = uc_x_range[2],
                            by = 0.02),
                    y = seq(from = uc_y_range[1],
                            to = uc_y_range[2],
                            by = 0.1))  # expand points to grid


# create an empty grid for whole of river
grd3 <- expand.grid(x = seq(from = all_x_range[1],
                            to = all_x_range[2],
                            by = 0.02),
                    y = seq(from = all_y_range[1],
                            to = all_y_range[2],
                            by = 0.1))  # expand points to grid


# create coord lists
b_list_low <- list(x = buff_spdf_low$x, y = buff_spdf_low$y)
b_list_up <- list(x = buff_spdf_up$x, y = buff_spdf_up$y)
b_list_all <- list(x = buff_spdf_all$x, y = buff_spdf_all$y)

# create irregular grid following buffered bottom - lower
C_grd_low <- grd1[!splancs::inout(grd1, b_list_low),]
coordinates(C_grd_low) <- ~x + y
gridded(C_grd_low) <- TRUE

# create irregular grid following buffered bottom - upper
C_grd_up <- grd2[!splancs::inout(grd2, b_list_up),]
coordinates(C_grd_up) <- ~x + y
gridded(C_grd_up) <- TRUE

# create irregular grid following buffered bottom - all river
C_grd_all <- grd3[!splancs::inout(grd3, b_list_all),]
coordinates(C_grd_all) <- ~x + y
gridded(C_grd_all) <- TRUE


## Save out sysdtat.rda
usethis::use_data(chl_brk, S_sitesdf, C_sitesdf, S_oxy_locs, C_oxy_locs,
                  S_bottom, S_bottom_nar, C_bottom_open, C_bottom_weir,
                  S_grd_all, S_grd_nar, C_grd_low, C_grd_up, C_grd_all,
                  reclass_matrices, internal = TRUE, overwrite = TRUE)
