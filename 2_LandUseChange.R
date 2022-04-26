## %######################################################%##
#                                                          #
####               Land Use Change                      ####
#                                                          #
## %######################################################%##

{
  require(terra)
  require(dplyr)
  require(flexsdm)
  require(here)
  require(ggplot2)
  require(viridis)
  require(raster)
  require(FedData)
}

## %######################################################%##
#                                                          #
####       Additional Data Layers for Cropping          ####
#                                                          #
## %######################################################%##

# California floristic provinces
cfp <- terra::vect("./1_Inputs/4_Shapefiles/JepsonRegions/JepsonRegions_CFP.shp")

# Environmental variables - Current conditions
env <- file.path('1_Inputs/2_Predictors/1_Current') %>%
  list.files(., full.names = T, pattern = '.tif$') %>%
  terra::rast()


## %######################################################%##
#                                                          #
####                      ICLUS Data                    ####
#                                                          #
## %######################################################%##

# Integrated Climate Land Use Scenario (ICLUS)
# https://www.epa.gov/iclus
# Developed by coupling human demographic growth model with a spatial allocation model (SERGoM) Theobald 2005
# Urban growth scenarios available are consistent with the IPCC Specieal Report on Emissions Scenarios (SRES) (Bierweagen et al. 2010)
# Decadal projections at 90-m spatial resolution (2000-2100)
# Comprehensive methods and metadata:
#   U.S. EPA. Updates To The Demographic And Spatial Allocation Models To Produce Integrated Climate And Land Use Scenarios
# (Iclus) (Final Report, Version 2). U.S. Environmental Protection Agency, Washington, DC, EPA/600/R-16/366F, 2017.
# Website download options
# Population projections (ICLUS v2.1.1)
# Land Use Projections under different climate change scenarios
# SSP2 RCP45 HadGEM2-ES & SSP2 RCP45 GISS-E2-R - seems like we could apply these to the RCP 45 climate projections
# SSP5 RCP85 HadGEM2-ES & SSP5 RCP85 GISS-E2-R - apply these to the RCP 85 projections

# Notes from Alex + Janet Meeting
# Reassign classes to: Developed, Agriculture, Other Land, Exurban development
# Proportional declines in suitability (this by nature will be somewhat arbitrary)
# I propose: Developed (100% reduction), Exurban development (50%), Agriculture (50%), Other Land (0%)

# Classes
# Natural water (Other Land)
# Reservoirs, canals (Other Land)
# Wetlands (Other Land)
# Recreation, conservation (Other Land)
# Timber (Other Land)
# Grazing (Other Land - this is basically shrubland)
# Pasture ()
# Cropland (Agriculture)
# Mining, barren land (Developed)
# Parks, golf courses (Developed)
# Exurban low density (Exurban)
# Exurban high density (Exurban)
# Suburban (Developed)
# Urban, low density (Developed)
# Urban, high density (Developed)
# Commercial (Developed)
# Industrial (Developed)
# Institutional (Other Land)
# Transportation (Developed)

# Folder with cropped land cover tiffs
dir <- list.dirs('./2_Outputs/8_Land_use/1_ICLUS/', recursive = FALSE)
names(dir) <- basename(dir)


f <- dir %>% list.files(full.names = TRUE,
                        recursive = TRUE,
                        pattern = '.tif$')

r <- terra::rast(f[1])


# Codes
codes <- seq(0,18,1)

# Class Names
class <-
  c('Natural water',
    'Reservoirs, canals',
    'Wetlands',
    'Recreations, conservation',
    'Timber',
    'Grazing',
    'Pasture',
    'Cropland',
    'Mining, barren land',
    'Parks, golf courses',
    'Exurban, low density',
    'Exurban, high density',
    'Suburban',
    'Urban, low density',
    'Urban, high density',
    'Commercial',
    'Industrial',
    'Institutional',
    'Transportation')

# Reclassifying where 1=natural, 0 = developed/agriculture, and .5 = exurban
recode <-
  c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, .5, .5, 0, 0, 0, 0, 0, 1, 0)

reclass <-
  c(
    'Natural',
    'Natural',
    'Natural',
    'Natural',
    'Natural',
    'Natural',
    'Developed',
    'Developed',
    'Developed',
    'Developed',
    'Exurban',
    'Exurban',
    'Developed',
    'Developed',
    'Developed',
    'Developed',
    'Developed',
    'Natural',
    'Developed'
  )

# raster attribute table
rat <- data.frame(ID = codes, class = class, recode = recode, reclass = reclass)


levels(r) <- rat
activeCat(r) <- 1

r1 <- r

par(mfrow = c(1, 2))
plot(r)
activeCat(r1) <- 3
plot(r1, col=c("darkgreen", "grey20", "yellow"))

# Reclassify
for(i in 1:length(dir)) {
  r <-
    paste0(dir[i], '/raw/') %>% list.files(pattern = '.tif$',
                                           full.names = T,
                                           recursive = T)
  
  r <- lapply(r, terra::rast)
  r <- lapply(r, subst, codes, recode)
  
  names(r) <- seq(2000, 2100, 10)
  
  nms <-
    paste0(names(r), ".tif") # extract raster names
  
  lgth <- length(nms)
  
  message("saving maps...")
  for (ii in 1:length(nms)) {
    terra::writeRaster(r[[ii]],
                       filename = file.path(dir[i], 'reclassified', nms[ii]),
                       overwrite = TRUE)
  }
}


############################################################
#                                                          #
#                Crop & Reproject LUC Data                 #
#                                                          #
############################################################

rcrs <-
  "G:/My Drive/Franklin_grant/project/data/landcover/ICLUS/ICLUS_v2_1_1_land_use_conus_ssp2_rcp45_hadgem2_es/ICLUS_v2_1_1_land_use_conus_2000.tif" %>%
  terra::rast() %>%
  terra::crs()

layer1 <- "1_Inputs/2_Predictors/1_Current/aet.tif" %>% terra::rast()
layer1[!is.na(layer1)] <- 1
layer2 <- terra::project(layer1, rcrs)
plot(layer1)
# plot(layer2)

res(layer2) / res(rast("2_Outputs/8_Land_use/1_ICLUS/ssp2_rcp45_hadgem2_es/reclassified/2000.tif"))


# List of directories
dirr <- "2_Outputs/8_Land_use/1_ICLUS/" %>%
  list.dirs() %>%
  grep("reclassified", ., value = TRUE)

for (d in 1:length(dirr)) {
  lfiles <- dirr[d] %>%
    list.files(pattern = ".tif$", full.names = T)
  luluc <- lfiles %>% lapply(terra::rast)
  luluc2 <- luluc %>% lapply(terra::aggregate, fact = 3, fun = mean) # aggregate with mean function
  luluc2 <- luluc2 %>% lapply(terra::project, terra::crs(layer1)) # project
  luluc2 <- luluc2 %>% lapply(crop, layer1) # crop
  luluc2 <- luluc2 %>% lapply(mask, mask = jep) # mask
  luluc2 <- luluc2 %>% lapply(terra::resample, layer1) # aggregate with resample() to have the exact resolution and extent
  luluc2 <- luluc2 %>% terra::rast()
  names(luluc2) <- gsub(".tif", "", basename(lfiles))
  # terra::compareGeom(luluc2[[1]], layer1)
  
  dirr[d] <- gsub("reclassified", "study_a", dirr[d])
  dir.create(dirr[d])
  
  luluc2[is.na(luluc2)] <- 0
  # plot(luluc2[[1:5]], col = rev(pals::ocean.haline(20))[-c(20:17)])
  writeRaster(luluc2, file.path(dirr[d], paste0(names(luluc2), ".asc")), overwrite=TRUE)
  
  # Interpolate for each period ??????????
  
  dirr[d] <- gsub("study_a", "interpolated", dirr[d])
  dir.create(dirr[d])
  years <- as.numeric(names(luluc2))
  for (i in 1:(length(years) - 1)) {
    interr <-
      flexsdm::interp(
        r1 = luluc2[[i]],
        r2 = luluc2[[i + 1]],
        y1 = years[i],
        y2 = years[i + 1]
      )
    writeRaster(interr, file.path(dirr[d], paste0(names(interr), ".asc")), overwrite=TRUE)
  }
}


# Repeat data for 2000 from 1995 to 1999
dirr <- "2_Outputs/8_Land_use/1_ICLUS/" %>%
  list.dirs() %>%
  grep("interpolated", ., value = TRUE)

dirr <- list.files(dirr, full.names = TRUE, pattern = "2000.asc$")

r <- dirr[1] %>% terra::rast()
r <- terra::rast(list(r, r, r, r, r))
names(r) <- 1995:1999
writeRaster(r, file.path(dirname(dirr[1]), paste0(names(r), ".asc")), overwrite=TRUE)

r <- dirr[2] %>% terra::rast()
r <- terra::rast(list(r, r, r, r, r))
names(r) <- 1995:1999
writeRaster(r, file.path(dirname(dirr[2]), paste0(names(r), ".asc")), overwrite=TRUE)



## %######################################################%##
#                                                          #
####               Prepare models with LULUC             ####
#                                                          #
## %######################################################%##

require(parallel)
require(doParallel)
require(raster)

lu_rcp_45 <- "./2_Outputs/8_Land_use/1_ICLUS/ssp2_rcp45_hadgem2_es/interpolated"  %>% list.files(pattern = ".asc$", full.names = TRUE)
lu_rcp_45 <- lu_rcp_45[1:91]
names(lu_rcp_45) <- gsub(".asc", "", basename(lu_rcp_45))

lu_rcp_85 <- "./2_Outputs/8_Land_use/1_ICLUS/ssp5_rcp85_hadgem2_es/interpolated" %>% list.files(pattern = ".asc$", full.names = TRUE)
lu_rcp_85 <- lu_rcp_85[1:91]
names(lu_rcp_85) <- gsub(".asc", "", basename(lu_rcp_85))


# Raw models
raw_models <- "./2_Outputs/9_Final_SDM/1_SDM_RAW" %>% list.files(full.names = TRUE)
names(raw_models) <- gsub(".zip$", "", basename(raw_models))

# Create directories
dir_save <- "./2_Outputs/9_Final_SDM/3_SDM_LULUC"

parallel::detectCores()
cl <- parallel::makeCluster(10, outfile="")
doParallel::registerDoParallel(cl)


# Process SDM
for (i in 1:length(raw_models)) {
  # message("Unzipping files", " sp num. ", i)
  # unzip(raw_models[i], exdir = dir_save)
  
  dir_save2 <- file.path(dir_save, names(raw_models[i]))
  dir_save2 <- dir_save2 %>% list.dirs(recursive = FALSE)
  files <- lapply(dir_save2, function(x) list.files(x, pattern = ".asc$", full.names = TRUE))
  names(files) <- basename(dir_save2)
  
  # Current condition
  lf <- files$`01_current`
  r <- raster::raster(lf)
  r <- r * raster::raster(lu_rcp_45[["1995"]])
  raster::writeRaster(r, lf, overwrite = TRUE, format = "ascii")
  
  # 02_cnrm_rcp45
  lf <- files$`02_cnrm_rcp45`
  message(dirname(lf[1]))
  r <- raster::raster(lf[1])
  r <- r * raster::raster(lu_rcp_45[["2055"]])
  raster::writeRaster(r, lf[1], overwrite = TRUE, format = "ascii")
  
  r <- raster::raster(lf[2])
  r <- r * raster::raster(lu_rcp_45[["2085"]])
  raster::writeRaster(r, lf[2], overwrite = TRUE, format = "ascii")
  
  
  # 04_hades_rcp45
  lf <- files$`04_hades_rcp45`
  message(dirname(lf[1]))
  r <- raster::raster(lf[1])
  r <- r * raster::raster(lu_rcp_45[["2055"]])
  raster::writeRaster(r, lf[1], overwrite = TRUE, format = "ascii")
  
  r <- raster::raster(lf[2])
  r <- r * raster::raster(lu_rcp_45[["2085"]])
  raster::writeRaster(r, lf[2], overwrite = TRUE, format = "ascii")
  
  
  # 03_cnrm_rcp85
  lf <- files$`03_cnrm_rcp85`
  message(dirname(lf[1]))
  r <- raster::raster(lf[1])
  r <- r * raster::raster(lu_rcp_85[["2055"]])
  raster::writeRaster(r, lf[1], overwrite = TRUE, format = "ascii")
  
  r <- raster::raster(lf[2])
  r <- r * raster::raster(lu_rcp_85[["2085"]])
  raster::writeRaster(r, lf[2], overwrite = TRUE, format = "ascii")
  
  
  # 05_hades_rcp85
  lf <- files$`05_hades_rcp85`
  message(dirname(lf[1]))
  r <- raster::raster(lf[1])
  r <- r * raster::raster(lu_rcp_85[["2055"]])
  raster::writeRaster(r, lf[1], overwrite = TRUE, format = "ascii")
  
  r <- raster::raster(lf[2])
  r <- r * raster::raster(lu_rcp_85[["2085"]])
  raster::writeRaster(r, lf[2], overwrite = TRUE, format = "ascii")
  
}
stopCluster(cl)



## %######################################################%##
#                                                           #
####      Prepare models with LULUC ONLY                 ####
#                                                           #
## %######################################################%##
require(parallel)
require(doParallel)
require(raster)

raw_models <- "./2_Outputs/9_Final_SDM/1_SDM_RAW" %>% list.files(full.names = TRUE)
names(raw_models) <- gsub(".zip$", "", basename(raw_models))
sp <- names(raw_models) %>% unique()

dir_sp <- file.path("./2_Outputs/9_Final_SDM/4_SDM_LULC_ONLY/", sp)
#sapply(dir_sp, dir.create)

folders <- c(
  "01_current", "02_hades_rcp45", "03_hades_rcp85"
)

#sapply(sapply(dir_sp, function(x) file.path(x, folders)), dir.create)

lu_rcp_45 <- "./2_Outputs/8_Land_use/1_ICLUS/ssp2_rcp45_hadgem2_es/interpolated"  %>% list.files(pattern = ".asc$", full.names = TRUE)
lu_rcp_45 <- lu_rcp_45[1:91]
names(lu_rcp_45) <- gsub(".asc", "", basename(lu_rcp_45))

lu_rcp_85 <- "./2_Outputs/8_Land_use/1_ICLUS/ssp5_rcp85_hadgem2_es/interpolated" %>% list.files(pattern = ".asc$", full.names = TRUE)
lu_rcp_85 <- lu_rcp_85[1:91]
names(lu_rcp_85) <- gsub(".asc", "", basename(lu_rcp_85))


# Create directories
dir_save <- "./2_Outputs/9_Final_SDM/4_SDM_LULC_ONLY"

# Raw models
raw_dir <- "./2_Outputs/9_Final_SDM/1_SDM_Raw"

n_list <-c(1:56,59:112,114:118)

# Process SDM
for (i in n_list) {
  
  raw_dir2 <- file.path(raw_dir, sp[i])
  raw_dir2 <- raw_dir2 %>% list.dirs(recursive = FALSE)
  files <- lapply(raw_dir2, function(x) list.files(x, pattern = ".asc$", full.names = TRUE))
  names(files) <- basename(raw_dir2)
  
  dir_save2 <- file.path(dir_save, sp[i])
  dir_save2 <- dir_save2 %>% list.dirs(recursive = FALSE)
  folders <- basename(dir_save2)
  
  # Current condition
  lf <- files$`01_current`
  r <- raster::raster(lf)
  r <- r * raster::raster(lu_rcp_45[["1995"]])
  
  terra::writeRaster(terra::rast(r),
                     paste0(dir_save2[1], '/', sp[i], '_1995.asc'),
                     overwrite = TRUE)
  

  # 02_hades_rcp45
  message(dir_save2[2])
  r <- raster::raster(lf)
  r <- r * raster::raster(lu_rcp_45[["2055"]])
  
  
  terra::writeRaster(terra::rast(r),
                     paste0(dir_save2[2], '/', sp[i], '_2055.asc'),
                     overwrite = TRUE)
  
  r <- raster::raster(lf)
  r <- r * raster::raster(lu_rcp_45[["2085"]])
  
  terra::writeRaster(terra::rast(r),
                     paste0(dir_save2[2], '/', sp[i], '_2085.asc'),
                     overwrite = TRUE)
  
  
  
  # 03_hades_rcp85
  message(dir_save2[3])
  r <- raster::raster(lf)
  r <- r * raster::raster(lu_rcp_85[["2055"]])
  
  terra::writeRaster(terra::rast(r),
                     paste0(dir_save2[3], '/', sp[i], '_2055.asc'),
                     overwrite = TRUE)
  
  r <- raster::raster(lf)
  r <- r * raster::raster(lu_rcp_85[["2085"]])
  
  
  terra::writeRaster(terra::rast(r),
                     paste0(dir_save2[3], '/', sp[i], '_2085.asc'),
                     overwrite = TRUE)
}
stopCluster(cl)

