##%######################################################%##
#                                                          #
####       Occurrence data cleaning and filtering       ####
#                                                          #
##%######################################################%##


# Packages and functions
{
  require(terra)
  require(dplyr)
  require(data.table)
  require(ggplot2)
  require(ggspatial)
  require(patchwork)
  require(tidyr)
  require(here)
}

# study area
studya <-  file.path("1_Inputs/2_Predictors/1_Current") %>%
  list.files(., full.names = T, pattern = 'aet') %>%
  terra::rast()
names(studya) <- 'study_area'
studya[studya > 0] <- 1

# CFP shapefile
cfp <- terra::vect("1_Inputs/4_Shapefiles/JepsonRegions/JepsonRegions_CFP.shp") %>%
  terra::project(crs(studya))

# Databases:
# raw occurrences database
db0 <-
  data.table::fread('1_Inputs/1_Occurrences/2_AllOccurrences/0_all_raw_data.gz') %>% tibble %>%
  relocate(date:elevation, .after = source) %>%
  relocate(county:locality, .after = elevation) %>%
  relocate(search_name:recorded_by, .after = elevation) %>%
  dplyr::select(
    -survey_type,
    -survey_date,
    -elevation,
    -sarea,
    -scientific_name,
    -latitude,
    -longitude
  ) %>%
  rename(species = search_name)

##%######################################################%##
#                                                          #
####                     Filtering                      ####
#                                                          #
##%######################################################%##

# Sort columns
ncell <-
  terra::cellFromXY(studya,
                    db0 %>% dplyr::select(longitude_m, latitude_m) %>% as.matrix)
db0 <- tibble(db0, ncell)

db0 <- db0 %>% dplyr::arrange(species, ncell, desc(year),  Grove_Name, data_base)

##%######################################################%##
#                                                          #
####                  Filter by cell                    ####
#                                                          #
##%######################################################%##
db0 <- db0 %>% group_by(species) %>% dplyr::filter(!duplicated(ncell))
db0 <- db0 %>% dplyr::select(-ncell)
db0 <- db0 %>% group_by()

##%######################################################%##
#                                                          #
####        Filter occurrence outside study area        ####
#                                                          #
##%######################################################%##

filt <- raster::extract(studya,
                        db0 %>% dplyr::select(longitude_m, latitude_m) %>% data.frame)
db0 <- db0 %>% dplyr::filter(!is.na(filt))

##%######################################################%##
#                                                          #
####             Filtering by institution               ####
#                                                          #
##%######################################################%##
# cal_inst <- data.table::fread("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/AuxiliaryData/CalInstitutions.txt") %>% tibble()
#
#
# ins_buf <- terra::vect(cal_inst, geom=c('x_m', 'y_m'), crs=crs(studya))
# ins_buf <- terra::buffer(ins_buf, 270*2)
# ins_buf$v <- 1
# ins_buf <- terra::aggregate(ins_buf, dissolve=TRUE, vars='v')
#
# filt <- terra::extract(ins_buf, terra::vect(db0, geom=c('longitude_m', 'latitude_m')))
# filt <- data.frame(filt)
# filt <- which(is.na(filt$id.y))
# dim(db0)
# db0 <- db0[filt,]
# dim(db0)

##%######################################################%##
#                                                          #
####                  Filter by year                    ####
#                                                          #
##%######################################################%##
table(is.na(db0$year))
dbNA <- db0 %>% dplyr::filter(is.na(year))
dbNA %>% count(species) %>% arrange(desc(n))
db0 <- db0 %>% dplyr::filter(year>=1980)
nrow(db0)


##%######################################################%##
#                                                          #
####              Filtering by occurrence               ####
####            geographical precision (for             ####
####         this filter NA will no be removed)         ####
#                                                          #
##%######################################################%##
db0 <- db0 %>% dplyr::filter(location_quality %in% c("medium", "high", ""))
dim(db0)


##%######################################################%##
#                                                          #
####       Filtering species occurrences based in       ####
####          CFP that CalFlora records falls           ####
#                                                          #
##%######################################################%##

db0 <- db0 %>% group_by()

ecofilt <-
  terra::extract(cfp[, 'Region'], db0 %>% terra::vect(., geom=c('longitude_m', 'latitude_m')))
db0 <- db0 %>%
  dplyr::mutate(ecofilt = ecofilt$Region) %>%
  tidyr::unite("sp_eco", species, ecofilt, remove = FALSE)

db0 %>% dplyr::filter(data_base == 'Calflora') %>% count(species) %>% arrange(n)

# Dataset with species and CFP regions
ecofilt <- db0 %>%
  dplyr::filter(data_base=='Calflora') %>%
  dplyr::distinct(species, ecofilt) %>%
  arrange(species) %>%
  tidyr::unite("sp_eco", 1:2, remove = FALSE)

db0 <- db0[db0$sp_eco%in%ecofilt$sp_eco, ]

db0 <- db0 %>% dplyr::select(-sp_eco)

db0[is.na(db0$ecofilt),]


##%######################################################%##
#                                                          #
####     detecting potential not wild occurrences       ####
#                                                          #
##%######################################################%##

non_wild <-
  c('botanic',
    'botanical',
    'zoo',
    'campus',
    'cultivated',
    'nursery',
    'garden',
    'campus')

filt <- grepl(paste(non_wild,collapse="|"),
              tolower(db0$locality))
table(filt) #38 occurrences have this words

db0 <- db0 %>% dplyr::mutate(wild=!filt)

# individualize each row
db0 <- db0 %>% group_by(data_base) %>% mutate(IDr=paste(data_base, 1:length(data_base)))
db0 <- db0 %>% group_by()


plot(cfp)
points(db0 %>% dplyr::filter(species=='Sequoia sempervirens') %>%
         dplyr::select(longitude_m, latitude_m), col='red', pch=19)


##%######################################################%##
#                                                          #
####                  Define absences                   ####
#                                                          #
##%######################################################%##

require(flexsdm)

# defining absence based on ecoregion (ecofilt)

cfp %>% plot
spp <- db0$species %>% unique


absences_eco <- db0 %>% dplyr::select(species, ecofilt) %>% distinct() %>% na.omit()
absences_nr <- db0 %>% count(species)

absences <- vector(mode = "list", length = length(spp))
names(absences) <- spp

for(i in 1:length(spp)) {
  print(i)
  f <-
    absences_eco %>% dplyr::filter(species == spp[i]) %>% pull(ecofilt)
  absences[[spp[i]]] <- db0 %>%
    dplyr::filter(species != spp[i]) %>%
    dplyr::filter(ecofilt %in% f) %>%
    dplyr::filter(survey == 'releve') %>%
    dplyr::select(-species) %>% unique()
}


absences <- dplyr::bind_rows(absences, .id='species')
absences <- absences %>% mutate(IDr=as.character(1:nrow(absences)))
absences$pr_ab <- 0


db2 <- db0 %>% dplyr::select(species, longitude_m, latitude_m, IDr, data_base, occ_id, survey, year, recorded_by, locality, ecofilt)
db2$pr_ab <- 1
db2 <- bind_rows(db2, absences)


##%######################################################%##
#                                                          #
####                  Filter by cell                    ####
#                                                          #
##%######################################################%##

# Sort columns
ncell <-
  terra::cellFromXY(studya,
                    db2 %>% dplyr::select(longitude_m, latitude_m) %>% as.matrix)
db2 <- tibble(db2, ncell)

db2 <- db2 %>% dplyr::arrange(species, desc(pr_ab), ncell, desc(year), data_base)


db2 <- db2 %>% group_by(species) %>% dplyr::filter(!duplicated(ncell))
db2 <- db2 %>% dplyr::select(-ncell)
db2 <- db2 %>% group_by()

##%######################################################%##
#                                                          #
####  Manually removing ecoregions for suspicious records  #
#                                                          #
##%######################################################%##

# Ecoregions: Central Western CA, Cascade Ranges, Great Valley, NorthWestern CA, Sierra Nevada, Southwestern CA

# Pinus torreyana only occurs in the Southwestern Ecoregion
# Sequoia sempervirens occurs Northwestern CA & Central Western CA
# Sequoiadendron giganteum occurs in the Sierra Nevada
# Arctostaphylos glandulosa does not occur in the Great Valley
# Arctostaphylos mewukka is only found in the Cascade Ranges and the Sierra Nevada
# Artemisia californica is found in the Southwestern CA and the Central Western CA
# Calocedrus decurens not found in Great Valley
# Calochortus albus not found in Great Valley
# Pinus coulteri only found in the Southwestern CA and the Central Western CA
# Pinus jeffreyi & Pinus ponderosa not found in Great Valley
# Prunus ilicifolia & virginiana is not found in Great Valley
# Quercus chrysolepis not found in Great Valley
# Quercus dumosa & engelmannii only found in Southwestern CA
# Quercus kelloggi not found in Great Valley
# Rhus integrifolia is only in Southwestern CA
# Salvia apiana not found in Great Valley
# Salvia leucophylla not found in Great Valley
# Salvia mellifera not found in Great Valley or Northwestern CA

db2 <- db2[!((
  db2$species == 'Pinus torreyana' &
    db2$ecofilt %in% c(
      'Central Western CA',
      'Cascade Ranges',
      'Great Valley',
      'NorthWestern CA',
      'Sierra Nevada'
    )
) |
  (
    db2$species == 'Sequoia sempervirens' &
      db2$ecofilt %in% c(
        'Cascade Ranges',
        'Great Valley',
        'Southwestern CA',
        'Sierra Nevada'
      )
  ) |
  (
    db2$species == 'Sequoiadendron giganteum' &
      db2$ecofilt %in% c(
        'Central Western CA',
        'Cascade Ranges',
        'Great Valley',
        'NorthWestern CA',
        'Southwestern CA'
      )
  )
|
  (
    db2$occ_id %in% c(
      'mg88502',
      'in:48341250',
      'in:52707945',
      'ce242',
      'mg106275',
      'in:66765840',
      'mg146395',
      'in:67487645',
      'in:35514109',
      'in:52063080',
      'in:52707945',
      'wb2234-435',
      'mg88502',
      'in:63803727',
      'mg76260',
      'io29244',
      'mg125195',
      'po65766',
      'gr2128',
      'bur6546',
      'ce870',
      'ce931',
      'svy1618',
      'svy1643',
      'ce592',
      'jgr23532',
      'ce943',
      'ce944',
      'svy1399',
      'svy1414',
      'svy1425',
      'svy1498',
      'svy1499',
      'svy1681',
      'bur6544',
      'in:3411859',
      'svy3030',
      'svy2979',
      'svy1452',
      'svy3563',
      'svy2897',
      'svy3182',
      'svy4075',
      'svy4111'
    )
  )
|
  (
    db2$species == 'Arctostaphylos glandulosa' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Arctostaphylos mewukka' &
      db2$ecofilt %in% c('Southwestern CA')
  )
|
  (
    db2$species == 'Artemisia californica' &
      db2$ecofilt %in% c('Sierra Nevada',
                         'Great Valley',
                         'NorthWestern CA')
  )
|
  (
    db2$species == 'Calocedrus decurrens' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Calochortus albus' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Pinus coulteri' &
      db2$ecofilt %in% c('Great Valley',
                         'NorthWestern CA')
  )
|
  (
    db2$species == 'Pinus jeffreyi' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Pinus ponderosa' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Prunus ilicifolia' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Prunus virginiana' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Quercus chrysolepis' &
      db2$ecofilt %in% c('Great Valley')
  )
| (
  db2$species == 'Quercus dumosa' &
    db2$ecofilt %in% c(
      'Central Western CA',
      'Cascade Ranges',
      'Great Valley',
      'NorthWestern CA',
      'Sierra Nevada'
    )
)
| (
  db2$species == 'Quercus engelmannii' &
    db2$ecofilt %in% c(
      'Central Western CA',
      'Cascade Ranges',
      'Great Valley',
      'NorthWestern CA',
      'Sierra Nevada'
    )
)
|
  (
    db2$species == 'Quercus kelloggii' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Rhus integrifolia' &
      db2$ecofilt %in% c('Great Valley', 'Central Western CA')
  )
|
  (
    db2$species == 'Salvia apiana' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Salvia leucophylla' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Salvia mellifera' &
      db2$ecofilt %in% c('Great Valley', 'NorthWestern CA')
  )
|
  (
    db2$species == 'Scutellaria californica' &
      db2$ecofilt %in% c('Great Valley')
  )
|
  (
    db2$species == 'Xylococcus bicolor' &
      db2$ecofilt %in% c('Central Western CA')
  )
), ]

plot(cfp)
points(db2 %>% dplyr::filter(species=='Sequoiadendron giganteum' & pr_ab ==1) %>%
         dplyr::select(longitude_m, latitude_m), col='red', pch=19)
points(db2 %>% dplyr::filter(species=='Sequoia sempervirens' & pr_ab ==1) %>%
         dplyr::select(longitude_m, latitude_m), col='blue', pch=19)
points(db2 %>% dplyr::filter(species=='Pinus torreyana' & pr_ab ==1) %>%
         dplyr::select(longitude_m, latitude_m), col='green', pch=19)



data.table::fwrite(db2, '1_Inputs/1_Occurrences/2_AllOccurrences/1_spp_pres_abs_cleaned.gz')

## %######################################################%##
#                                                          #
####             Create directory structure             ####
#                                                          #
## %######################################################%##

scenarios <- c(
  'cnrm_rcp45_2040_2069',
  'cnrm_rcp45_2070_2099',
  'cnrm_rcp85_2040_2069',
  'cnrm_rcp85_2070_2099',
  'hades_rcp45_2040_2069',
  'hades_rcp45_2070_2099',
  'hades_rcp85_2040_2069',
  'hades_rcp85_2070_2099'
)

# set up sdm directory
# dir <- sdm_directory(
#   main_dir = "directory path",
#   calibration_area = TRUE,
#   algorithm = c('glm', 'gam', 'gbm', 'raf', 'net', 'svm', 'gau', 'max'),
#   ensemble =  c('meansup', 'meanw', 'meanthr', 'mean', 'median'),
#   projections = scenarios,
#   threshold = TRUE,
#   return_vector = TRUE
# )
# 
# dir %>% head()

## %######################################################%##
#                                                          #
####                  Calibration area                  ####
#                                                          #
## %######################################################%##

# Occ database
occ <- data.table::fread('1_Inputs/1_Occurrences/2_AllOccurrences/1_spp_pres_abs_cleaned.gz') %>% tibble()
occ <- occ %>% dplyr::select(species, id = IDr, x = longitude_m, y = latitude_m, year, pr_ab)

# Jepson ecoregion
cfp <- terra::vect("1_Inputs/4_Shapefiles/JepsonRegions/JepsonRegions_CFP.shp")
cfp$Region = NULL

# Process species
sp <- unfilt_occ$species %>%
  unique() %>%
  sort()

for (i in 1:length(sp)) {
  x2 <-
    flexsdm::calib_area(
      data = occ[occ$species == sp[i], ],
      x = "x",
      y = "y",
      method = c("mask", cfp, "RegionCode")
    )
  x2$Region <- NULL
  terra::writeVector(x2, file.path(grep("Calibration", dir, value = TRUE), paste0(sp[i], ".shp")), overwrite = TRUE)
}

## %######################################################%##
#                                                          #
####                 Data partitioning                  ####
#                                                          #
## %######################################################%##
# Raterize CFP database
pred <-  file.path('1_Inputs/2_Predictors/1_Current') %>%
  list.files(., full.names = T, pattern = '.tif$') %>%
  terra::rast()
pred <- homogenize_na(pred)

cfp_r <- rasterize(cfp, pred$aet, field = "RegionCode")
cfp_r[] <- as.numeric(cfp_r[])

# Extract ecoregion value
occ$region <- terra::extract(cfp_r, occ[, c("x", "y")])[, "RegionCode"]
occ <- occ %>% filter(!is.na(region))

# Block partition for species with > 30 presences
nr <- occ %>% dplyr::filter(pr_ab == 1) %>% count(species)
sp <- nr$species[nr$n > 30]

part <- list()
for (i in 1:length(sp)) {
  print(i)
  part[[i]] <-
    flexsdm::part_sblock(
      env_layer = pred[[!names(pred) %in% "terrain"]],
      data = occ %>% filter(species == sp[i]),
      x = "x",
      y = "y",
      pr_ab = "pr_ab",
      n_part = 4, # four blocks
      min_res_mult = 50,
      max_res_mult = 300,
      num_grids = 30,
      min_occ = 15,
      prop = 0.1
    )
}
names(part) <- sp

# Species without good partition
sp <- names(which(sapply(part, function(x) is.null(names(x)))))

# Three folds
for (i in 1:length(sp)) {
  print(i)
  part[[sp[i]]] <-
    flexsdm::part_sblock(
      env_layer = pred[[!names(pred) %in% "terrain"]],
      data = occ %>% filter(species == sp[i]),
      x = "x",
      y = "y",
      pr_ab = "pr_ab",
      n_part = 3, # three blocks
      min_res_mult = 50,
      max_res_mult = 300,
      num_grids = 30,
      min_occ = 15,
      prop = 0.5
    )
}


# Species without good partition
sp_n <- names(which(sapply(part, function(x) is.null(names(x)))))

part <- part[sapply(part, length) == 3]

# saver block partition details and rasters
best_part_info <- bind_rows(lapply(part, function(x) x$best_part_info), .id = "species")
occ_part <- bind_rows(lapply(part, function(x) x$part), .id = "species")
occ_part %>%
  group_by(species, .part) %>%
  count() %>%
  filter(n == min(n))
data.table::fwrite(best_part_info, file.path("1_Inputs/1_Occurrences/3_OccPartition/best_part_info_block.gz"))
data.table::fwrite(occ_part, file.path("1_Inputs/1_Occurrences/3_OccPartition/occ_part_block.gz"))
blocks <- lapply(part, function(x) {
  get_block(pred$aet, x$grid)
})
blocks <- rast(blocks)
names(blocks) <- names(part)
writeRaster(blocks, file.path("1_Inputs/1_Occurrences/3_OccPartition/1_SpatialBlocks/", paste0('part_',names(blocks), ".tif")))


#### Partition for those species with nr$n>20 <= 30

# "Poa stebbinsii", "Ribes lasianthum", "Arctostaphylos rainbowensis", "Artemisia rothrockii", "Pinus torreyana"

# Band partition
sp <- nr$species[nr$n > 20 & nr$n <= 30]
sp <- c(sp_n, sp)
part_ba <- as.list(rep(NA, length(sp)))
names(part_ba) <- sp
# pb <- progress_bar$new(total = length(sp))
# plot_res(pred[[1]], 300)
for (i in 1:length(sp)) {
  print(i)
  # pb$tick()
  part_ba[[sp[i]]] <-
    part_sband(
      env_layer = pred[[!names(pred) %in% "terrain"]],
      data = occ %>% filter(species == sp[i]),
      x = "x",
      y = "y",
      pr_ab = "pr_ab",
      type = "lat",
      n_part = 2, # Two bands
      min_bands = 3,
      max_bands = 30,
      min_occ = 10,
      prop = 0.5
    )
}
names(part_ba) <- sp

sp_n <- names(which(sapply(part_ba, function(x) is.null(names(x)))))
# "Arctostaphylos rainbowensis" - one of Santiago's species

# save band partition details and rasters
part_ba <- part_ba[sapply(part_ba, length) == 3]
best_part_info_ba <- bind_rows(lapply(part_ba, function(x) x$best_part_info), .id = "species")
occ_part_ba <- bind_rows(lapply(part_ba, function(x) x$part), .id = "species")
occ_part_ba %>%
  group_by(species, .part) %>%
  count() %>%
  filter(n == min(n))
data.table::fwrite(best_part_info_ba, file.path("1_Inputs/1_Occurrences/3_OccPartition/best_part_info_band.gz"))
data.table::fwrite(occ_part_ba, file.path("1_Inputs/1_Occurrences/3_OccPartition/occ_part_band.gz"))
blocks <- lapply(part_ba, function(x) {
  get_block(pred$aet, x$grid)
})
blocks <- rast(blocks)
names(blocks) <- names(part_ba)
writeRaster(blocks, file.path("1_Inputs/1_Occurrences/3_OccPartition/1_SpatialBlocks/", paste0('part_',names(blocks), ".tif")))


best_part_info <- data.table::fread(file.path("1_Inputs/1_Occurrences/3_OccPartition/best_part_info_block.gz"))


## %######################################################%##
#                                                          #
####                     4-Modeling                     ####
#                                                          #
## %######################################################%##

{
  require(dplyr)
  require(terra)
  require(flexsdm)
  require(here)
  require(progress)
  require(raster)
  require(ggplot2)
  require(kernlab)
  require(foreach)
}



# wmean ensemble
wmean_ens <- function(m, w, thr){
  m <- terra::weighted.mean(m, w)
  m_thr <- m*(m>=thr)
  m <- c(m, m_thr)
  names(m) <- c('meanw', 'max_sens_spec')
  return(m)
}

# mean ensemble
mean_ens <- function(m, thr){
  m <- terra::mean(m)
  m_thr <- m*(m>=thr)
  m <- c(m, m_thr)
  names(m) <- c('mean', 'max_sens_spec')
  return(m)
}

# median ensemble
median_ens <- function(m, thr){
  m <- terra::median(m)
  m_thr <- m*(m>=thr)
  m <- c(m, m_thr)
  names(m) <- c('median', 'max_sens_spec')
  return(m)
}

# mean of superior models
meansup_ens <- function(m, w, thr) {
  m <- terra::mean(m[[which(w >= mean(w))]])
  m_thr <- m*(m>=thr)
  m <- c(m, m_thr)
  names(m) <- c('meansup', 'max_sens_spec')
  return(m)
}

# mean threshold ensemble
meanthr_ens <- function(m, thr_m, thr){
  for(r in 1:terra::nlyr(m)){
    m[[r]][m[[r]]<thr_m[r]] <- 0
  }
  m <- terra::mean(m)
  m_thr <- m*(m>=thr)
  m <- c(m, m_thr)
  names(m) <- c('meanthr', 'max_sens_spec')
  return(m)
}

# mean thereshold ensemble for only one raster
meanthr_ens2 <- function(m, thr_m, thr){
  m[m<thr_m] <- 0
  m <- terra::mean(m)
  m_thr <- m*(m>=thr)
  m <- c(m, m_thr)
  names(m) <- c('meanthr', 'max_sens_spec')
  return(m)
}


##%######################################################%##
#                                                          #
####             Prepare data for modeling              ####
####           species with > 50 occurrences            ####
#                                                          #
##%######################################################%##
# Read occurrences databases
occ <- data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_block.gz') %>% tibble()
occ2 <-  data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_band.gz') %>% tibble()

occ <- bind_rows(occ, occ2)

# Environmental variables - Current conditions
env <- file.path('1_Inputs/2_Predictors/1_Current') %>%
  list.files(., full.names = T, pattern = '.tif$') %>%
  terra::rast()
env <- homogenize_na(env)
env_names <- names(env)


env_fut <- file.path("1_Inputs/2_Predictors/2_Projection/") %>% list.dirs(recursive = F)
names(env_fut) <- env_fut
factor <- 'terrain'
env_fut <- as.list(env_fut)

for(i in 1:length(env_fut)) {
  env_fut[[i]] <- env_fut[[i]] %>%
    list.files(., pattern = ".tif$", full.names = TRUE) %>%
    terra::rast()
}
names(env_fut) <- basename(names(env_fut))

# Extract environmental variables
occ <- sdm_extract(occ, x = "x", y = "y", env_layer = env, filter_na = TRUE)

# Count number of presences
n_occ <- occ %>% dplyr::filter(pr_ab==1) %>% pull(species) %>% table() %>% sort()
sp <- names(n_occ)[n_occ>50]

##%######################################################%##
#                                                          #
####             Prepare  data for modeling             ####
####           species with <=50 occurrences            ####
#                                                          #
##%######################################################%##
# Read occurrences databases
occ <- data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_block.gz') %>% tibble()
occ2 <-  data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_band.gz') %>% tibble()

occ <- bind_rows(occ, occ2)

# Environmental variables - Current conditions
# topo <-
#   "G:/My Drive/Franklin_grant/project/data/Topo_hetero/topo_hetero.tif" %>% terra::rast() %>%
#   project(env[[1]]) %>%
#   terra::crop(cfp) %>%
#   terra::mask(cfp)
# 
# writeRaster(topo, "G:/My Drive/Franklin_grant/project/data/Topo_hetero/topo_hetero_cfp.tif")

topo <- terra::rast("G:/My Drive/Franklin_grant/project/data/Topo_hetero/topo_hetero_cfp.tif")

env <- file.path('1_Inputs/2_Predictors/1_Current') %>%
  list.files(., full.names = T, pattern = '.tif$') %>%
  terra::rast()
env <- homogenize_na(env)
env$terrain <- NULL
env <- rast(list(env, topo))
plot(env)
env_names <- names(env)

# Future predictors
env_fut <- "./1_Inputs/2_Predictors/2_Projection" %>% list.dirs(recursive = F)
names(env_fut) <- env_fut
env_fut <- as.list(env_fut)

for(i in 1:length(env_fut)) {
  env_fut[[i]] <- env_fut[[i]] %>%
    list.files(., pattern = ".tif$", full.names = TRUE) %>%
    terra::rast()
  env_fut[[i]]$terrain <- NULL
  env_fut[[i]] <- rast(list(env_fut[[i]], topo))
}

names(env_fut) <- basename(names(env_fut))


# Extract env conditions
occ <- sdm_extract(occ, x = "x", y = "y", env_layer = env, filter_na = TRUE)

# Count number of presences
n_occ <- occ %>% dplyr::filter(pr_ab==1) %>% pull(species) %>% table() %>% sort()
sp <- names(n_occ)[n_occ<=50]


##%######################################################%##
#                                                          #
####                 Loop for modeling                  ####
#                                                          #
##%######################################################%##
modeled <- "2_Outputs/1_Current/Ensemble/meanthr/1_con" %>% list.files() %>% gsub(".tif$", "",.)

perf_dir <- here("2_Outputs/0_Model_performance/")

# 
n_list <- c(91, 78, 84, 60, 47)

for (i in n_list) {
  message(paste("Modeling sp", i, sp[i]))
  pa <- occ %>% dplyr::filter(species == sp[i])
  
  #### Boosted regression trees ####
  m_gbm <- flexsdm::tune_gbm(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    grid = expand.grid(
      n.trees = seq(10, 200, 30),
      shrinkage = seq(0.1, 2, 0.5),
      n.minobsinnode = seq(1, 15 , 3)
    ),
    thr = "max_sens_spec",
    metric = "TSS",
    n_cores = 4 # length(pa$.part %>% unique())
  )
  
  if(length(m_gbm)>1){
    h <- m_gbm$hyper_performance
    
    p1 <- h %>% ggplot(aes(x = n.trees, y = AUC_mean,
                           col = as.factor(n.minobsinnode)), group = shrinkage) +
      geom_line() +
      facet_wrap(. ~ shrinkage) +
      theme(legend.position = "bottom") +
      labs(x = "n.minobsinnode")
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_gbm_auc',  '.png')), dpi=200)
    
    p1 <- h %>% ggplot(aes(x = n.trees, y = TSS_mean,
                           col = as.factor(n.minobsinnode)), group = shrinkage) +
      geom_line() +
      facet_wrap(. ~ shrinkage) +
      theme(legend.position = "bottom") +
      labs(x = "n.minobsinnode")
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_gbm_tss',  '.png')), dpi=200)
    
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_gbm.txt")))
  }
  
  #### Neural Network ####
  m_net <- tune_net(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    grid = expand.grid(
      size = (2:length(env_names)),
      decay = c(seq(0.01, 1, 0.05), 1, 3, 4, 5, 10)
    ),
    thr = "max_sens_spec",
    metric = "TSS",
    n_cores = 4
  )
  
  if(length(m_net)>1){
    h <- m_net$hyper_performance
    
    p1 <- h %>% ggplot(aes(x = decay, y = AUC_mean, col = as.factor(size))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_net_auc',  '.png')), dpi=200)
    
    p1 <- h %>% ggplot(aes(x = decay, y = TSS_mean, col = as.factor(size))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_net_tss',  '.png')), dpi=200)
    
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_net.txt")))
  }
  
  
  #### Random forest ####
  m_raf <- tune_raf(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    grid = expand.grid(mtry = seq(1, length(env_names), 1)),
    thr = "max_sens_spec",
    metric = "TSS",
    n_cores = 4
  )
  
  if(length(m_raf)>1){
    h <- m_raf$hyper_performance
    
    h %>% ggplot(aes(x = mtry, y = AUC_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_raf_auc',  '.png')), dpi=200)
    h %>% ggplot(aes(x = mtry, y = TSS_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_raf_tss',  '.png')), dpi=200)
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_raf.txt")))
  }
  
  
  #### Support Vector Machine ####
  m_svm <- tune_svm2(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    grid = expand.grid(C = seq(2, 80, 20),
                       sigma = seq(0.001, 0.2, 0.05))
    ,
    thr = "max_sens_spec",
    metric = "TSS",
    n_cores = 4
  )
  
  if(length(m_svm)>1){
    h <- m_svm$hyper_performance
    
    h %>% ggplot(aes(x = sigma, y = AUC_mean, col = as.factor(C))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_svm_auc',  '.png')), dpi=200)
    h %>% ggplot(aes(x = sigma, y = TSS_mean, col = as.factor(C))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i],' hyp_svm_tss',  '.png')), dpi=200)
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_svm.txt")))
  }
  
  #### Generalized Additive Model ####
  n_t <- flexsdm:::n_training(data = pa, partition = ".part")
  
  candidate_k <- 20
  while(any(n_t < flexsdm:::n_coefficients(data = pa, predictors = env_names[!env_names %in% factor], predictors_f = factor, k = candidate_k))){
    candidate_k <- candidate_k-3
  }
  
  # candidate_k <- 20
  # while(any(n_t < flexsdm:::n_coefficients(data = pa, predictors = env_names[!env_names %in% factor], k = candidate_k))){
  #   candidate_k <- candidate_k-3
  # }
  
  m_gam <- fit_gam(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    thr = "max_sens_spec",
    k = candidate_k
  )
  
  #### Generalized Linear Models ####
  if(sum(pa$pr_ab==1)>=100){
    m_glm <- fit_glm(
      data = pa,
      response = "pr_ab",
      # predictors = env_names[c(1:8, 10)],
      predictors = env_names[!env_names %in% factor],
      predictors_f = factor,
      partition = ".part",
      thr = "max_sens_spec",
      poly = 2
    )
  }
  
  
  # compile models objects
  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]
  
  ##### Ensemble ####
  m_ensemble <-
    flexsdm::fit_ensemble(
      lapply(models, get),
      ens_method = c("meanw", "meanthr", "mean", "meansup", "median"),
      thr_model = "max_sens_spec",
      metric = "TSS",
      thr = "max_sens_spec"
    )
  
  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]
  
  # Model performance
  performance <- flexsdm::sdm_summarize(lapply(models, function(x) {
    if (length(get(x)) > 0) {
      get(x)
    }
  }))
  
  readr::write_tsv(x = performance, file = here(perf_dir, paste0(sp[i], "_models_performance.txt")))
  
  
  ##%######################################################%##
  #                                                          #
  ####             Predict individual models              ####
  #                                                          #
  ##%######################################################%##
  models <- models[models!="m_ensemble"]
  
  message("Predicting models for species ", sp[i], " ", i)
  
  models_object <- lapply(models, function(x) {
    get(x)
  })
  
  prd <-
    sdm_predict(
      models = models_object,
      pred = env,
      thr = c("max_sens_spec"),
      con_thr = TRUE,
      clamp = TRUE,
      pred_type = 'cloglog',
      predict_area = NULL
    )
  
  for (mm in 1:length(prd)) {
    terra::writeRaster(
      prd[[mm]],
      file.path(
        './2_Outputs/1_Current',
        'Algorithm',
        names(prd[mm]),
        '1_con',
        paste0(sp[i], '.tif')
      )
      ,
      overwrite = TRUE
    )
  }
  
  
  ##### Model projection #####
  for (f in 1:length(env_fut)) {
    print(f)
    prd <-
      sdm_predict(
        models = models_object,
        pred = env_fut[[f]],
        thr = c("max_sens_spec"),
        con_thr = TRUE,
        clamp = TRUE,
        pred_type = 'cloglog',
        predict_area = NULL
      )
    
    for (mm in 1:length(prd)) {
      terra::writeRaster(
        prd[[mm]],
        file.path(
          './2_Outputs/2_Projection',
          names(env_fut[f]),
          'Algorithm',
          names(prd[mm]),
          '1_con',
          paste0(sp[i], '.tif')
        )
        ,
        overwrite = TRUE
      )
    }
    rm(prd)
  }
  
  rm(list = grep("m_", ls(), value = TRUE))
}


##%######################################################%##
#                                                          #
####                 Predict ensemble                   ####
#                                                          #
##%######################################################%##

# Read occurrences databases
occ <-
  data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_block.gz') %>% tibble()
occ2 <-
  data.table::fread('1_Inputs/1_Occurrences/3_OccPartition/occ_part_band.gz') %>% tibble()
occ <- bind_rows(occ, occ2)

n_occ <- occ %>% dplyr::filter(pr_ab==1) %>% pull(species) %>% table() %>% sort()
sp <- names(n_occ)[n_occ>50]

perf_dir <- "./2_Outputs/0_Model_performance/"

env_fut <- "./1_Inputs/2_Predictors/2_Projection" %>% list.dirs(recursive = F)
names(env_fut) <- env_fut
env_fut <- as.list(env_fut)
names(env_fut) <- basename(names(env_fut))


for (i in n_list) {
  message(paste("Predicting ensemble models for sp", i, sp[i]))
  pa <- occ %>% dplyr::filter(species == sp[i])
  
  ens_perf <-
    paste0(perf_dir, sp[i], "_models_performance.txt") %>%
    readr::read_tsv(col_types = readr::cols())
  
  #models_perf <- gsub("m_", "", ens_perf) %>% sort()
  
  
  filt_perf <- ens_perf$AUC_mean >= 0.7 &
    ens_perf$thr_value != 0 &
    ens_perf$thr_value != 1
  
  models_perf <- ens_perf$model[filt_perf]
  models_perf <-
    models_perf[!models_perf %in% c("meanthr", "meanw", "median", "meansup", "mean")]
  
  thr <-
    ens_perf$thr_value[ens_perf$model %in% c("meanthr", "meanw", "median", "meansup", "mean")]
  
  w <- ens_perf %>% dplyr::filter(model %in% models_perf) %>% arrange(model) %>% pull(TSS_mean)
  
  thr_mod <- ens_perf$thr_value[ens_perf$model %in% models_perf]
  
  dd <- './2_Outputs/1_Current/Algorithm' %>%
    list.dirs(., recursive = FALSE)
  dd <- dd[basename(dd) %in% models_perf]
  dd <-
    lapply(dd, function(x)
      list.files(
        x,
        pattern =  sp[i],
        recursive = TRUE,
        full.names = TRUE
      )) %>% unlist() %>%
    terra::rast()
  dd <- dd[[grep("max_sens_spec", names(dd), invert = TRUE)]]
  dd <- dd[[names(dd) %>% sort]]
  
  # Predict ensemble
  prd <- list(meanthr_ens(m = dd, thr_m = thr_mod, thr = thr[1]),
              wmean_ens(m = dd, w = w, thr = thr[2]),
              median_ens(m = dd, thr = thr[3]),
              meansup_ens(m = dd, w = w, thr = thr[4]),
              mean_ens(m = dd, thr = thr[5]))
  names(prd) <- c("meanthr", "meanw", "median", "meansup", "mean")
  
  for (mm in 1:length(prd)) {
    terra::writeRaster(
      prd[[mm]],
      file.path(
        './2_Outputs/1_Current/Ensemble',
        names(prd[mm]),
        '1_con',
        paste0(sp[i], '.tif')
      )
      ,
      overwrite = TRUE
    )
  }
  
  rm(prd)
  
  
  ##### Model projection #####
  for (f in 1:length(env_fut)) {
    print(f)
    
    dd <-
      file.path('./2_Outputs/2_Projection',
                names(env_fut[f]),
                'Algorithm') %>%
      list.dirs(., recursive = FALSE)
    dd <- dd[basename(dd) %in% models_perf]
    dd <-
      lapply(dd, function(x)
        list.files(
          x,
          pattern =  sp[i],
          recursive = TRUE,
          full.names = TRUE
        )) %>% unlist() %>%
      terra::rast()
    dd <- dd[[grep("max_sens_spec", names(dd), invert = TRUE)]]
    dd <- dd[[names(dd) %>% sort]]
    
    # for Rhus integrifolia remove net model (very low suitability values)
    dd <- dd[[c('glm', 'svm')]]
    
    # Predict ensemble
    prd <- list(meanthr_ens(m = dd, thr_m = thr_mod, thr = thr[1]),
                wmean_ens(m = dd, w = w, thr = thr[2]),
                median_ens(m = dd, thr = thr[3]),
                meansup_ens(m = dd, w = w, thr = thr[4]),
                mean_ens(m = dd, thr = thr[5]))
    names(prd) <- c("meanthr", "meanw", "median", "meansup", "mean")
    
    for (mm in 1:length(prd)) {
      terra::writeRaster(
        prd[[mm]],
        file.path(
          './2_Outputs/2_Projection/',
          names(env_fut[f]),
          'Ensemble',
          names(prd[mm]),
          '1_con',
          paste0(sp[i], '.tif')
        )
        ,
        overwrite = TRUE
      )
    }
    
    rm(prd)
    
  }
}


##%######################################################%##
#                                                          #
####      Crop models and correct overprediction        ####
#                                                          #
##%######################################################%##

# Occupied habitat (before correcting for land use patterns)

# Directories for saving outputs
dir_save <- "./2_Outputs/9_Final_SDM/1_SDM_Raw/" %>% list.dirs(., full.names = TRUE, recursive = FALSE)
dir_save <- file.path(dir_save, "01_current")
names(dir_save) <- basename(dirname(dir_save))

dir_sp <-  "./2_Outputs/9_Final_SDM/1_SDM_Raw/" %>% list.dirs(., full.names = TRUE, recursive = FALSE)
names(dir_sp) <- basename(dir_sp)

# California floristic provinces
cfp <- terra::vect("./1_Inputs/4_Shapefiles/JepsonRegions/JepsonRegions_CFP.shp")

# Read occurrences databases and extract ecoregions
occ <- data.table::fread(
  file.path(
    getwd() ,
    "1_Inputs/1_Occurrences/2_AllOccurrences/1_spp_pres_abs_cleaned.gz"
  )
) %>%
  tibble() %>%
  dplyr::rename(x = longitude_m, y = latitude_m) %>%
  dplyr::select(species, x, y, pr_ab)

# p <-
#   occ %>% dplyr::filter(pr_ab == 1) %>% terra::vect(., geom = c("x", "y"), crs = terra::crs(cfp))
#
# eco <- terra::extract(cfp, p)
# eco <- tibble(eco[-1])
# p$RegionCode <- eco$RegionCode
# rm(eco)

# terra::writeVector(
#   p,
#   "1_Inputs/1_Occurrences/2_AllOccurrences/2_spp_pres_abs_cleaned_w_eco.shp",
#   overwrite = TRUE
# )


p <- terra::vect("1_Inputs/1_Occurrences/2_AllOccurrences/2_spp_pres_abs_cleaned_w_eco.shp")


# List of model - current conditions
list_f <- "./2_Outputs/1_Current/Ensemble/meanthr/1_con" %>% list.files(full.names = TRUE)
names(list_f) <- basename(list_f) %>% gsub(".tif$", "", .)

# Crop models and constrain predictions
for(i in 1:length(list_f)){
  
  print(names(list_f[i]))
  
  # Remove suitability values outside calibration area
  r <- list_f[i] %>% terra::rast()
  r <- r[[2]] # Thresholded model
  p_sp <- p[p$species == names(list_f[i])]
  eco <- p_sp$RegionCode %>% unique()
  eco <- c('NW', 'CaR', 'GV', 'SN', 'CW') # manual for blue oak
  cfp2 <- cfp[cfp$RegionCode %in% eco]
  r1 <- terra::mask(r, cfp2)
  r1[is.na(r1) & !is.na(r)] <- 0 # replace NA with 0
  r1[is.na(r1)] <- 0
  
  # Correct overprediction based on pres approach
  occ_models_sp <-
    occ %>% dplyr::filter(species == names(list_f[i])) %>% na.omit()
  r2 <- flexsdm::msdm_posteriori(
    records = occ_models_sp,
    x = "x",
    y = "y",
    pr_ab = "pr_ab",
    cont_suit = r1,
    method = "pres", # some species this did not work (Scutellaria californica and Erigeron petrophilus)
   # method = "mcp", # less strict
    thr = "max_sens_spec"
  )
  r2 <- r2[[1]] #select the continuous model
  r2[is.na(r2)] <- 0 # replace NA with 0
  
  png(file.path(dir_sp[names(list_f[i])], paste0("sp_current_range.png")), width=10, height=10, units="in", res=300)
  par(mfrow=c(1, 2))
  plot(r1, main='Original habitat', cex=1.3)
  plot(vect(occ_models_sp %>% filter(pr_ab ==1), geom = c('x', 'y')), add = TRUE, alpha = .5, cex = .2)
  plot(cfp, add = TRUE)
  plot(r2, main='Occupied habitat', cex=1.3)
  plot(cfp, add = TRUE)
  dev.off()
  
  # Save model in 01_current  folder
  terra::writeRaster(
    r2,
    filename = file.path(dir_save[names(list_f[i])], paste0(names(list_f[i]), "_1995.asc")),
    overwrite = TRUE
  )
}


