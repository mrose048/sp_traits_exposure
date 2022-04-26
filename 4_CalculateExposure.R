##%######################################################%##
#                                                          #
####                Calculate exposure                  ####
#                                                          #
##%######################################################%##

{
  require(terra)
  require(dplyr)
  require(flexsdm)
  require(here)
  require(ggplot2)
  require(viridis)
  require(tidyr)
  require(landscapemetrics)
  require(raster)
  require(landscapetools)
  require(data.table)
  require(parallel)
  require(doParallel)
  require(raster)
}


# California floristic provinces
cfp <- terra::vect("./1_Inputs/4_Shapefiles/JepsonRegions/JepsonRegions_CFP.shp")


##%######################################################%##
#                                                          #
####        Calculating spp exposure                    ####
####                                                    ####
#                                                          #
##%######################################################%##

# function for calculating exposure
exp_fun <-
  function(m1, m2, thr, res = 270, sp, scen, rcp, model, year) {
    names(m1) <- 'current'
    
    names(m2) <- 'future' # full dispersal
    
    m2_mask <-
      terra::mask(m2, m1, maskvalue = c(0, NA), updatevalue = 0) # null dispersal future
    
    m1[is.na(m1)] <-
      0 # replace NA with 0
    m2[is.na(m2)] <-
      0 # replace NA with 0
    
    nd_count <-
      terra::global(m2_mask, sum)[[1]] # null dispersal suitability-weighted cell #
    c_count <-
      terra::global(m1, sum)[[1]] # current suitability_weighted cell #
    
    exp_df <-
      tibble(
        sp,
        model,
        year,
        scen,
        rcp,
        current_cells = c_count,
        nd_future_cells = nd_count,
        current_km2 = c_count * (res * res) / 1000000,
        nd_future_km2 = nd_count * (res * res) / 1000000,
        nd_change_prop = (nd_count-c_count)/c_count,
        exposure = ((nd_count-c_count)/c_count) * -1
      )
    
    return(exp_df)
  }

thr <-
  readr::read_tsv("./2_Outputs/0_Model_performance/00_model_performance.txt") %>%
  filter(species != 'Mimulus cardinalis')

thr2 <- readr::read_tsv("./2_Outputs/0_Model_performance/00_model_performance_pres_only.txt")

thr <- bind_rows(thr, thr2)

# species to process
lu_dir <- "./2_Outputs/9_Final_SDM/2_SDM_LU" %>% list.dirs(., recursive = FALSE)
sapply(lu_dir, function(x) length(list.files(x))) %>% unique() # ok
names(lu_dir) <- basename(lu_dir)
sp <- names(lu_dir)

# Climate change only
for (i in 1:length(sp)) {
  
  print(sp[i])
  
  dir_save2 <- file.path("./2_Outputs/9_Final_SDM/2_SDM_LU", sp[i])
  dir_save2 <- dir_save2 %>% list.dirs(recursive = FALSE)
  files <- lapply(dir_save2, function(x) list.files(x, pattern = ".asc$", full.names = TRUE))
  names(files) <- basename(dir_save2)
  
  # threshold for species meanthr model
  thr2 <- thr %>%
    dplyr::filter(model == 'meanthr' ,
                  species == sp[i]) %>%
    dplyr::select(thr_value)
  
  
  # Current condition
  lf <- files$`01_current`
  r <- terra::rast(lf)
  
  # 02_cnrm_rcp45
  lf <- files$`02_cnrm_rcp45`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  cnrm45 <- rbind(
    exposure(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2055'
    ),
    exposure(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  
  # 04_hades_rcp45
  lf <- files$`04_hades_rcp45`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  hades45 <- rbind(
    exposure(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2055'
    ),
    exposure(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  # 03_cnrm_rcp85
  lf <- files$`03_cnrm_rcp85`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  cnrm85 <- rbind(
    exposure(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2055'
    ),
    exposure(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  # 05_hades_rcp85
  lf <- files$`05_hades_rcp85`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  hades85 <- rbind(
    exposure(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2055'
    ),
    exposure(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  
  exp_df <- bind_rows(cnrm45, hades45, cnrm85, hades85)
  df <- exp_df %>%
    dplyr::mutate(change_driver = 'climate_change')
  
  readr::write_tsv(df,
                   file =
                     file.path(
                       "./2_Outputs/3_Exposure/data",
                       paste0(sp[i], "_sp_and_exp_traits_cc.txt")
                     ))
  
}

# Land use change only
for (i in 1:length(sp)) {
  
  print(sp[i])
  
  dir_save2 <- file.path("./2_Outputs/9_Final_SDM/4_SDM_LULC_ONLY", sp[i])
  dir_save2 <- dir_save2 %>% list.dirs(recursive = FALSE)
  files <- lapply(dir_save2, function(x) list.files(x, pattern = ".asc$", full.names = TRUE))
  names(files) <- basename(dir_save2)
  
  # threshold for species meanthr model
  thr2 <- thr %>%
    dplyr::filter(model == 'meanthr' ,
                  species == sp[i]) %>%
    dplyr::select(thr_value)
  
  
  # Current condition
  lf <- files$`01_current`
  r <- terra::rast(lf)
  
  c_range_traits <- range_traits(m1 = r, thr = thr2$thr_value, sp = sp[i], model = 'meanthr')
  
  
  # 02_hades_rcp45
  lf <- files$`02_hades_rcp45`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  crs(r55) <- crs(r)
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  crs(r85) <- crs(r)
  
  hades45 <- rbind(
    exp_fun(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2055'
    ),
    exp_fun(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  
  
  # 03_hades_rcp85
  lf <- files$`03_hades_rcp85`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  crs(r55) <- crs(r)
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  crs(r85) <- crs(r)
  
  hades85 <- rbind(
    exp_fun(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2055'
    ),
    exp_fun(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  
  exp_df <- bind_rows(hades45, hades85)
  df <- left_join(c_range_traits, exp_df, by = c('sp', 'model')) %>%
    dplyr::mutate(change_driver = 'lulc')
  
  readr::write_tsv(df,
                   file =
                     file.path(
                       "./2_Outputs/3_Exposure/data",
                       paste0(sp[i], "_sp_and_exp_traits_lulc_only.txt")
                     ))
  
}

# Climate + land use change 
for (i in 1:length(sp)) {
  
  print(sp[i])
  
  dir_save2 <- file.path("./2_Outputs/9_Final_SDM/3_SDM_LULUC", sp[i])
  dir_save2 <- dir_save2 %>% list.dirs(recursive = FALSE)
  files <- lapply(dir_save2, function(x) list.files(x, pattern = ".asc$", full.names = TRUE))
  names(files) <- basename(dir_save2)
  
  # threshold for species meanthr model
  thr2 <- thr %>%
    dplyr::filter(model == 'meanthr' ,
                  species == sp[i]) %>%
    dplyr::select(thr_value)
  
  
  # Current condition
  lf <- files$`01_current`
  r <- terra::rast(lf)
  
  # 02_cnrm_rcp45
  lf <- files$`02_cnrm_rcp45`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  cnrm45 <- rbind(
    exp_fun(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2055'
    ),
    exp_fun(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  
  # 04_hades_rcp45
  lf <- files$`04_hades_rcp45`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  hades45 <- rbind(
    exp_fun(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2055'
    ),
    exp_fun(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp45',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  # 03_cnrm_rcp85
  lf <- files$`03_cnrm_rcp85`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  cnrm85 <- rbind(
    exp_fun(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2055'
    ),
    exp_fun(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'cnrm',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  # 05_hades_rcp85
  lf <- files$`05_hades_rcp85`
  message(dirname(lf[1]))
  
  r55 <- terra::rast(grep('2055', lf, value = TRUE))
  r85 <- terra::rast(grep('2085', lf, value = TRUE))
  
  hades85 <- rbind(
    exp_fun(
      m1 = r,
      m2 = r55,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2055'
    ),
    exp_fun(
      m1 = r,
      m2 = r85,
      thr = thr2$thr_value,
      sp = sp[i],
      scen = 'hades',
      rcp = 'rcp85',
      model = 'meanthr',
      year = '2085'
    )
  )
  
  
  exp_df <- bind_rows(cnrm45, hades45, cnrm85, hades85)
  df <- left_join(c_range_traits, exp_df, by = c('sp', 'model')) %>%
    dplyr::mutate(change_driver = 'climate_change + lulc')
  
  readr::write_tsv(df,
                   file =
                     file.path(
                       "./2_Outputs/3_Exposure/data",
                       paste0(sp[i], "_sp_and_exp_traits_cc_lulc.txt")
                     ))
  
}


##%######################################################%##
#                                                          #
####            Exposure Data                           ####
#                                                          #
##%######################################################%##

exp <- "./2_Outputs/3_Exposure/data/" %>% list.files(pattern = ".txt", full.names = TRUE)
exp <-
  lapply(exp,
         readr::read_tsv,
         col_types = list(
           .default = "d",
           sp = "c",
           model = 'c',
           scen = 'c',
           rcp = 'c',
           change_driver = 'c'
         ))
exp <- bind_rows(exp)

exp_df <- exp %>%
  mutate(range_km2 = c_ca/100) %>%
  dplyr::filter(model == 'meanthr') 

# checking number of observations by species
exp_df$sp %>% table  # all okay

# write merged exposure data to file
data.table::fwrite('2_Outputs/3_Exposure/exposure_data_raw.csv')




